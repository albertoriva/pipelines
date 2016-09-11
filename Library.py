# (c) 2015, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import os
import time
import math

import Utils
from Lines import Line
from SampleCollection import SampleCollection
from dmaputils import initializeContrasts, writeBSMAPscripts, writeMCOMPqsub, diffmethSummary, drawBEDplots, drawDMCplots, getCscallReport
import Table

### Testing

class testLine(Line):
    name = "Test"

    def Report(self):
        ACT = self.actor
        out = ACT.out
        ACT.scene(1, "Testing")
        Table.fileToTable("/ufrc/riva/ariva/tests/data.csv", "t1", stream=out, headerFromFile=True, visibleRows=10, maxrows=20)
        return True
        
### Samples manager

class Samples(Line):
    """This object handles initializing the SampleCollection, checking that all files exist, fixing paths."""
    name = "SamplesManager"

    def Setup(self):
        ACT = self.actor
        ACT.sc = SampleCollection(ACT.Conf)
        SC = ACT.sc

        ## Fix paths to fastq files in each readset and store base name.
        for rs in SC.readsets:
            rs['left'] = ACT.checkPath(rs['left'])
            rs['lbase'] = ACT.setFileExt(os.path.split(rs['left'])[1], "", remove=[".fastq", ".gz"]) # basename of left-side reads
            if rs['paired']:
                rs['right'] = ACT.checkPath(rs['right'])
                rs['rbase'] = ACT.setFileExt(os.path.split(rs['right'])[1], "", remove=[".fastq", ".gz"]) # basename of right-side reads
        return True

    def Verify(self):
        return self.actor.sc.verify()

    def Report(self):
        ACT = self.actor
        SC = ACT.sc

        ACT.scene("Input data")
        ACT.reportf("The following table summarizes the samples, conditions, and contrasts in this analysis. A readset is either a single fastq file or a pair of fastq files (for paired-end sequencing).")
        tbl1 = Table.ScrollingTable(id='tbl1', caption='Summary of input data')
        tbl1.startHead()
        tbl1.addHeaderRow(["Category", "Data"])
        tbl1.startBody()
        tbl1.addSectionRow("Summary of input data")
        tbl1.addRowHeader("Contrasts:")
        tbl1.addRow([ ", ".join([ contr['test'] + " vs. " + contr['control'] for contr in SC.contrasts])])
        tbl1.addRowHeader("Conditions:")
        tbl1.addRow([ ", ".join([ cond['name'] for cond in SC.conditions])])
        for cond in SC.conditions:
            condsamples = SC.conditionSamples(cond)
            tbl1.addRowHeader(cond['name'], rowspan=len(condsamples))
            for smp in condsamples:
                tbl1.addRow([ "{} ({} readsets)".format(smp['name'], len(smp['readsets'])) ])
        tbl1.toHTML(ACT.out)
        return True


### Samples manager (methylation)

class MethSamples(Samples):
    """This object handles initializing the SampleCollection for the dmap pipeline."""
    name = "SamplesManager (methylation)"

    def PreExecute(self):
        ACT = self.actor
        initializeContrasts(ACT)
        return True

### Samples manager (rnaseq)

class RNAseqSamples(Samples):
    name = "SamplesManager (RNAseq)"

    def PreExecute(self):
        ACT = self.actor
        SC = ACT.sc

        for smp in SC.samples:
            smp['mix'] = ACT.getConf('mix', section=smp['name'], default='1')

        return True
    
### Sequence trimming with Trimmomatic or Sickle (with FastQC before and after)

class Trimmer(Line):
    """This objects perform sequence trimming, using either trimmomatic or sickle. Use
the 'prog' property to select between the two ('trimmomatic' or 'sickle'). This also
performs quality control with FastQC before and after the trimming. Use the 'qcbefore'
and 'qcafter' properties to specify the output directories for the respective steps, or
set either one to False to disable that QC step)."""
    name = "Trimming and QC"
    trim = "trimmomatic"
    qcbefore = "QC-pretrim/"
    qcafter = "QC-trim/"

    def Setup(self):
        self.trim = Utils.dget('prog', self.properties, self.trim)
        self.qcbefore = Utils.dget('qcbefore', self.properties, self.qcbefore)
        self.qcafter = Utils.dget('qcafter', self.properties, self.qcafter)
        return True

    def Execute(self):
        nqc = 0
        ntrim = 0
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log
        props = self.properties
        if self.qcbefore:
            ACT.mkdir(self.qcbefore)
        if self.qcafter:
            ACT.mkdir(self.qcafter)

        for rs in SC.readsets:
            if rs['paired']:
                rs['leftpretrim'] = rs['left']
                rs['rightpretrim'] = rs['right']
                rs["left"]   = rs['lbase'] + ".trim.paired.fastq.gz"
                rs["right"]   = rs['rbase'] + ".trim.paired.fastq.gz"
                rs["trunpair1"] = rs['lbase'] + ".trim.unpaired.fastq.gz"
                rs["trunpair2"] = rs['rbase'] + ".trim.unpaired.fastq.gz"

                if self.qcbefore:
                    LOG.log("Performing QC on `{}' and `{}'", rs['leftpretrim'], rs['rightpretrim'])
                    rs['fastqcbeforeleft'] = Utils.fastqcPath(self.qcbefore, rs['leftpretrim'])
                    rs['fastqcbeforeright'] = Utils.fastqcPath(self.qcbefore, rs['rightpretrim'])
                    if not self.dry:
                        ACT.submit("fastqc.qsub {} {}".format(rs['leftpretrim'], self.qcbefore), done="fqc.@.done")
                        ACT.submit("fastqc.qsub {} {}".format(rs['rightpretrim'], self.qcbefore), done="fqc.@.done")
                        nqc += 2

                if self.trim == 'trimmomatic':
                    LOG.log("Calling trimmomatic on `{}' and `{}'", rs['leftpretrim'], rs['rightpretrim'])
                    LOG.log("Trimmomatic main output to `{}' and `{}'", rs['left'], rs['right'])
                    if not self.dry:
                        tr = ACT.submit("trimmomatic.qsub {} {} {} {} {} {}".format(
                                rs["leftpretrim"], rs["rightpretrim"],
                                rs["left"], rs["trunpair1"], 
                                rs["right"], rs["trunpair2"]),
                                        done="trim.@.done")
                        ntrim += 1
                elif props['prog'] == 'sickle':
                    print "*** not implemented yet"

                if self.qcafter:
                    LOG.log("Performing QC on `{}' and `{}'", rs['left'], rs['right'])
                    rs['fastqcafterleft'] = Utils.fastqcPath(self.qcafter, rs['left'])
                    rs['fastqcafterright'] = Utils.fastqcPath(self.qcafter, rs['right'])
                    if not self.dry:
                        ACT.submit("fastqc.qsub {} {}".format(rs["left"], self.qcafter), after=tr, done="fqc.@.done")
                        ACT.submit("fastqc.qsub {} {}".format(rs["right"], self.qcafter), after=tr, done="fqc.@.done")
                        nqc += 2
            else:
                rs['leftpretrim'] = rs['left']
                rs["left"]   = rs['lbase'] + ".trim.paired.fastq.gz"
                rs["trunpair1"] = rs['lbase'] + ".trim.unpaired.fastq.gz"

                if self.qcbefore:
                    LOG.log("Performing QC on `{}'", rs['leftpretrim'])
                    rs['fastqcbeforeleft'] = Utils.fastqcPath(self.qcbefore, rs['leftpretrim'])
                    if not self.dry:
                        ACT.submit("fastqc.qsub {} {}".format(rs['leftpretrim'], self.qcbefore), done="fqc.@.done")
                        nqc += 1

                if self.trim == 'trimmomatic':
                    LOG.log("Calling trimmomatic on `{}'", rs['leftpretrim'])
                    LOG.log("Trimmomatic main output to `{}'", rs['left'])
                    if not self.dry:
                        tr = ACT.submit("trimmomatic-se.qsub {} {}".format(rs["leftpretrim"], rs["left"]),
                                        done="trim.@.done")
                        ntrim += 1
                elif props['prog'] == 'sickle':
                    print "*** not implemented yet"

                if self.qcafter:
                    LOG.log("Performing QC on `{}'", rs['left'])
                    rs['fastqcafterleft'] = Utils.fastqcPath(self.qcafter, rs['left'])
                    if not self.dry:
                        ACT.submit("fastqc.qsub {} {}".format(self.qcafter, rs["left"]), after=tr, done="fqc.@.done")
                        nqc += 1
        LOG.log("Waiting for {} qc jobs and {} trimmer jobs", nqc, ntrim)
        return ACT.wait([("fqc.@.done", nqc), ("trim.@.done", ntrim)])

    def Report(self):
        ACT = self.actor
        SC = ACT.sc

        ACT.scene("Trimming and quality control")
        ACT.reportf("""The input sequences were trimmed using <b>{}</b>. Quality control was performed before and after trimming using <b>FastQC</b>. The following table provides links to the 
quality control reports before and after trimming, as well as the number of reads in the trimmed files.""".format(self.trim))
        prereadcounts = Utils.fileToDict(ACT.fastqCountsPreTrim)
        readcounts = Utils.fileToDict(ACT.fastqCounts)
        tbl2 = Table.ScrollingTable(id='tbl2', align="LRLRL", caption='Number of reads in input files and links to QC reports.')
        tbl2.startHead()
        tbl2.addHeaderRow(["Sample", "Readset", "Reads before trim", "QC before trim", "Reads after trim", "QC after trim"])
        tbl2.startBody()
        for smp in SC.samples:
            tbl2.addRowHeader(smp['name'], rowspan=len(smp['readsets']) + 1)
            for rs in smp['readsets']:
                if rs['paired']:
                    tbl2.addRow([ rs['name'], 
                                  Utils.fmt(prereadcounts[rs['name']]),
                                  rs['fastqcbeforeleft'] + "<br>" + rs['fastqcbeforeright'],
                                  Utils.fmt(readcounts[rs['name']]),
                                  rs['fastqcafterleft'] + "<br>" + rs['fastqcafterright']
                                  ])
                else:
                    tbl2.addRow([ rs['name'], 
                                  Utils.fmt(prereadcounts[rs['name']]),
                                  rs['fastqcbeforeleft'],
                                  Utils.fmt(readcounts[rs['name']]),
                                  rs['fastqcafterleft']
                                  ])
            tbl2.addRow([ "<b>" + smp['name'] + "</b> (total)", Utils.fmt(prereadcounts[smp['name']]),"&nbsp;", Utils.fmt(readcounts[smp['name']]), 
                          Utils.pct(readcounts[smp['name']], prereadcounts[smp['name']])])
        tbl2.toHTML(ACT.out)
        return True

### STAR aligner

class StarAligner(Line):
    """This object implements the STAR aligner."""
    name = "STAR Aligner"
    version = ""

    def Setup(self):
        ACT = self.actor
        SC = ACT.sc
        
        ## Get name of STAR index for each sample
        for s in SC.samples:
            s['staridx'] = ACT.getConf("staridx", s['name']) or ACT.staridx
        return True

    def Verify(self):
        ACT = self.actor
        version = ACT.shell("module load star; STAR --version")
        if version.startswith("STAR"):
            self.version = version
            return True
        else:
            return False

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nstar = 0
        for smp in SC.samples:
            # Create an output directory for each readset                                                                                                                        
            for rs in smp['readsets']:
                outdir = rs["name"]
                ACT.mkdir(outdir)
                rs["bam"] = rs["name"] + "/Aligned.sortedByCoord.out.bam" # Name of BAM file produced by STAR                                                               
                if not self.dry:
                    if rs['paired']:
                        LOG.log("Running STAR aligner on `{}' and `{}'", rs["left"], rs["right"])
                        ACT.submit("star.qsub {} {} {} {}".format(smp["staridx"], rs["left"], rs["right"], outdir), done="star.@.done")
                    else:
                        LOG.log("Running STAR aligner on `{}'", rs["left"])
                        ACT.submit("star.qsub {} {} {}".format(smp["staridx"], rs["left"], outdir), done="star.@.done")
                    nstar += 1
        return ACT.wait(("star.@.done", nstar))

### STAR aligner for transcriptome

class StarAlignerTx(Line):
    """This object implements the STAR aligner."""
    name = "STAR Aligner (transcriptome)"
    version = ""

    def Setup(self):
        ACT = self.actor
        SC = ACT.sc
        
        ## Get name of STAR index for each sample
        for s in SC.samples:
            s['staridx'] = ACT.getConf("staridx", s['name']) or ACT.staridx
        return True

    def Verify(self):
        ACT = self.actor
        version = ACT.shell("module load star; STAR --version")
        if version.startswith("STAR"):
            self.version = version
            return True
        else:
            return False

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nstar = 0
        for smp in SC.samples:
            # Create an output directory for each sample
            outdir = smp["name"]
            ACT.mkdir(outdir)
            smp["bam"] = smp["name"] + "/Aligned.toTranscriptome.out.bam" # Name of transcriptome BAM file produced by STAR
            smp["genomebam"] = smp["name"] + "/Aligned.sortedByCoord.out.bam" # Name of genome BAM file produced by STAR
            lfq = []
            rfq = []
            for rs in smp['readsets']:
                lfq.append(rs['left'])
                if rs['paired']:
                    rfq.append(rs['right'])
            lfq = ",".join(lfq)
            rfq = ",".join(rfq)
            if rfq == '':
                LOG.log("Running STAR aligner on `{}'", lfq)
            else:
                LOG.log("Running STAR aligner on `{}' and `{}'", lfq, rfq)
            if not self.dry:
                ACT.submit("star-exp.qsub {} {} {} {}".format(smp["staridx"], lfq, rfq, outdir), done="star.@.done")
                nstar += 1
        return ACT.wait(("star.@.done", nstar))

    def Report(self):
        ACT = self.actor
        SC = ACT.sc

        readcounts = Utils.fileToDict(ACT.fastqCounts)

        ACT.scene("Alignment to transcriptome")
        ACT.reportf("""The input sequences were aligned to the transcriptome using <b>{}</b>. The following table reports the number of
alignments to the genome and the transcriptome for each sample. Please note that the number of alignments will in general be higher than the number of
reads because the same read may align to multiple isoforms of the same gene. The WIG files can be uploaded to the <A href='http://genome.ucsc.edu/'>UCSC
Genome Browser</A> as <A href='http://genome.ucsc.edu/cgi-bin/hgCustom'>custom tracks</A>.""", self.version)
        geAlignments = Utils.fileToDict(ACT.genomeCounts)
        txAlignments = Utils.fileToDict(ACT.afterCounts)
        tbl3 = Table.ScrollingTable(id='tbl3', align="RRRRC", caption='Number of alignments to genome and transcriptome.')
        tbl3.startHead()
        tbl3.addHeaderRow(["Sample", "Genome alignments", "Genome alignment rate", "Transcriptome alignments", "Transcriptome alignment rate", "WIG file"])
        tbl3.startBody()
        for smp in SC.samples:
            sn = smp['name']
            tbl3.addRowHeader(sn)
            tbl3.addRow([ Utils.fmt(Utils.dget(sn, geAlignments, default=0)), 
                          Utils.f2dd(1.0 * Utils.dget(sn, geAlignments, default=0) / readcounts[sn]),
                          Utils.fmt(txAlignments[sn]),
                          Utils.f2dd(1.0 * txAlignments[sn] / readcounts[sn]),
                          Utils.linkify(smp['wigfile']) ])
        tbl3.toHTML(ACT.out)
        return True

### Fastq Counter

class FASTQcounter(Line):
    """Count the number of reads in the FASTQ files for all readsets. All FASTQ counting processes run in parallel,
and results are collected in the file indicated by the `outfile' property. The number of reads in each readset
is also stored in the `propname' field of the readset, and the total number of reads in each sample is stored in the
`propname' field of the sample. If `delay' is True, processing is delayed until the end of the pipeline run."""
    name = "FASTQ counter"
    nfq = 0
    prefix = ""
    outfile = None
    propname = None
    delay = None
    
    def Setup(self):
        self.outfile = Utils.dget('outfile', self.properties, self.outfile)
        self.propname = Utils.dget('propname', self.properties, self.propname)
        self.delay = Utils.dget('delay', self.properties, self.outfile)
        return True

    def writeQsub(self, qsubname, fastqfilenames, outprefix):
        with open(qsubname, "w") as out:
            out.write("""#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load dibig_tools

fastqfile={}

L=$1
fq=`lines $fastqfile $L`
countseqs.py $fq > {}-$L.out
""".format(fastqfilenames, outprefix))

    def combineCounts(self, SC, prefix):
        idx = 1
        good = True
        out = open(self.outfile, "w") if self.outfile else None
        # self.actor.shell("ls -l {}*".format(prefix))
        try:
            for smp in SC.samples:
                smptotal = 0
                for rs in smp['readsets']:
                    cntfile = "{}-{}.out".format(prefix, idx)
                    line = Utils.safeReadLineFromFile(cntfile, delay=2, maxtries=30)
                    if line:
                        pc = line.split("\t")
                        c = int(pc[1])
                        if c == 0:
                            self.status = "Warning: empty FASTQ file {}".format(cntfile)
                            good = False
                        else:
                            smptotal += c
                            if self.propname:
                                rs[self.propname] = c
                            if out:
                                out.write("{}\t{}\n".format(rs['name'], c))
                    else:
                        self.status = "Warning: empty FASTQ file {}".format(cntfile)
                        good = False
                    idx += 1
                if self.key:
                    smp[self.propname] = smptotal
                if out:
                    out.write("{}\t{}\n".format(smp['name'], smptotal))
        finally:
            if out:
                out.close()
        return good

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        id = Utils.id_generator(10)
        self.prefix = id
        fastqnames = id + "-names.txt"
        qsub = id + ".qsub"
        self.done = id + ".@.done"

        # Store the names of all FASTQ files in a file
        with open(fastqnames, "w") as out:
            for smp in SC.samples:
                for rs in smp['readsets']:
                    out.write("{}\n".format(rs['left']))
                    self.nfq += 1
        
        LOG.log("Submitting array job to compute number of reads in {} FASTQ files.", self.nfq)
        # Write the qsub script
        self.writeQsub(qsub, fastqnames, id)
        if not self.dry:
            for i in range(1, self.nfq + 1):
                ACT.submit("{} {}".format(qsub, i), done=self.done)
            if not self.delay:
                ACT.wait((self.done, self.nfq))
                if self.outfile:
                    LOG.log("Collecting FASTQ counts into {}.", self.outfile)
                if self.propname:
                    LOG.log("Storing FASTQ counts in key {}.", self.propname)
                return self.combineCounts(SC, id) # Signal error if one or more FASTQ files are empty
        return True

    def PostExecute(self):
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log

        good = True
        if not self.dry and self.delay:
            ACT.wait((self.done, self.nfq))
            if self.outfile:
                LOG.log("Collecting FASTQ counts into {}.", self.outfile)
            if self.propname:
                LOG.log("Storing FASTQ counts in key {}.", self.propname)
            good = self.combineCounts(SC, self.prefix) # Signal error if one or more FASTQ files are empty

        LOG.log("Deleting temporary count files.")
        ACT.shell("rm -f {}*".format(self.prefix))
        return True

### BAM Counter

class BAMcounter(Line):
    """Count the number of alignments in the BAM files for all samples. All BAM counting processes run in parallel,
and results are collected in the file indicated by the `outfile' property. The number of alignments in each readset
is also stored in the `key' field of the sample. If `delay' is True, processing is delayed until the end of the
pipeline run."""
    name = "BAM counter"
    nbams = 0
    prefix = ""
    source = "bam"              # property of sample containing name of BAM file
    outfile = None
    propname = None
    delay = None
    
    def Setup(self):
        self.outfile = Utils.dget('outfile', self.properties, self.outfile)
        self.propname = Utils.dget('propname', self.properties, self.propname)
        self.delay = Utils.dget('delay', self.properties, self.outfile)
        self.source = Utils.dget('source', self.properties, self.source)
        return True

    def writeQsub(self, qsubname, bamfilenames, outprefix):
        with open(qsubname, "w") as out:
            out.write("""#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load bamtools

bamfile={}

L=$1
bam=`lines $bamfile $L`
bamtools count -in $bam > {}-$L.out
""".format(bamfilenames, outprefix))

    def combineCounts(self, SC, prefix):
        idx = 1
        good = True
        out = open(self.outfile, "w") if self.outfile else None
        # self.actor.shell("ls -l {}*".format(prefix))
        try:
            for smp in SC.samples:
                cntfile = "{}-{}.out".format(prefix, idx)
                c = Utils.safeReadIntFromFile(cntfile, default=0, delay=2, maxtries=30)
                if c == 0:
                    self.status = "Warning: empty BAM file count in {}".format(cntfile)
                    good = False
                else:
                    if self.key:
                        smp[self.key] = c
                    if out:
                        out.write("{}\t{}\n".format(smp['name'], c))
                    idx += 1
        finally:
            if out:
                out.close()
        return good

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        id = Utils.id_generator(10)
        self.prefix = id
        bamnames = id + "-names.txt"
        qsub = id + ".qsub"
        self.done = id + ".@.done"

        # Store the names of all BAM files in a file
        with open(bamnames, "w") as out:
            for smp in SC.samples:
                out.write("{}\n".format(smp[self.source]))
                self.nbams += 1
        
        LOG.log("Submitting array job to compute number of reads in {} BAM files.", self.nbams)
        # Write the qsub array script
        self.writeQsub(qsub, bamnames, id)
        if not self.dry:
            for i in range(1, self.nbams+1):
                ACT.submit("{} {}".format(qsub, i), done=self.done)
            if not self.delay:
                ACT.wait((self.done, self.nbams))
                if self.outfile:
                    LOG.log("Collecting BAM counts into {}.", self.outfile)
                if self.propname:
                    LOG.log("Storing BAM counts in key {}.", self.propname)
                return self.combineCounts(SC, id) # Signal error if one or more BAM files are empty
        return True

    def PostExecute(self):
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log

        good = True
        if not self.dry and self.delay:
            ACT.wait((self.done, self.nbams))
            if self.outfile:
                LOG.log("Collecting BAM counts into {}.", self.outfile)
            if self.propname:
                LOG.log("Storing BAM counts in key {}.", self.propname)
            good = self.combineCounts(SC, self.prefix) # Signal error if one or more BAM files are empty
        
        LOG.log("Deleting temporary count files.")
        ACT.shell("rm -f {}*".format(self.prefix))
        return True

### Remove duplicates

class Markdup(Line):
    """Use Picard to remove duplicates and add the RG tag to the BAM file."""
    name = "Mark duplicates"
    
    def picardMarkDups(self, infile, outfile, metrics=None):
        cmdline = """picard.qsub MarkDuplicates \
I={} \
O={} \
CREATE_INDEX=true \
REMOVE_DUPLICATES=TRUE \
""".format(infile, outfile)
        if metrics:
            cmdline += """METRICS_FILE={} \
""".format(metrics)
        return cmdline

    def picardAddReplGroups(self, infile, outfile, rglb, rgsm):
        return """picard.qsub AddOrReplaceReadGroups \
INPUT={} \
OUTPUT={} \
CREATE_INDEX=TRUE \
RGPU=CARA \
RGID=CARA \
RGLB={} \
RGPL=illumina \
RGSM={} \
""".format(infile, outfile, rglb, rgsm)

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc
        
        npicard = 0
        for smp in SC.samples:
            smp['origbam'] = smp['bam']
            smp['markdup'] = smp['name'] + ".markdup.bam"
            smp['bam'] =   smp['name'] + ".addrg.bam"
            LOG.log("Removing duplicates and adding RG: {} -> {}", smp['origbam'], smp['bam'])
            if not self.dry:
                pic1 = ACT.submit(self.picardMarkDups(smp['origbam'], smp['markdup'], metrics=smp['name'] + "-metrics.txt"))
                ACT.submit(self.picardAddReplGroups(smp['markdup'], smp['bam'], smp['name'], smp['name']), after=pic1, done="picard.@.done")
                npicard += 1
        return ACT.wait(("picard.@.done", npicard))

    def PostExecute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        LOG.log("Deleting unnecessary BAM files.")
        if False: # not self.dry:
            for smp in SC.samples:
                Utils.safeRemoveFile(smp['origbam'])
                Utils.safeRemoveFile(smp['markdup'])
        return True

### BAM merger 

class BAMmerger(Line):
    """Merge sets of BAM files with bamtools."""

    name = "BAM merger"
    byCondition = False
    removeOriginal = False

    def Setup(self):
        self.byCondition = Utils.dget('byCondition', self.properties, self.byCondition)
        self.removeOriginal = Utils.dget('remove', self.properties, self.removeOriginal)
        return True

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nmerge = 0
        for smp in SC.samples:
            smp['bam'] = smp['name'] + ".bam"
            LOG.log("Merging BAM files for sample {} into {}", smp['name'], smp['bam'])
            cmdline = "generic.qsub bamtools merge -out " + smp['bam']
            for rs in smp['readsets']:
                cmdline += " -in " + rs['bam']
            cmdline += " module:bamtools"
            if not self.dry:
                ACT.submit(cmdline, done='merge.@.done')
                nmerge += 1
        ACT.wait(('merge.@.done', nmerge))

        if self.byCondition:
            nmerge = 0
            for cond in SC.conditions:
                cond['bam'] = cond['name'] + ".bam"
                LOG.log("Merging BAM files for condition {} into {}", cond['name'], cond['bam'])
                cmdline = "generic.qsub bamtools merge -out " + cond['bam']
                for sb in SC.conditionBAMs(cond['name']):
                    cmdline += " -in " + sb
                cmdline += " module:bamtools"
                if not self.dry:
                    ACT.submit(cmdline, done='cmerge.@.done')
                    nmerge += 1
            ACT.wait(('cmerge.@.done', nmerge))

        return True

### BAM merger 

class BAMconcatenator(Line):
    """Concatenate sets of BAM files with samtools."""

    name = "BAM concatenator"
    byCondition = False
    removeOriginal = False

    def Setup(self):
        self.byCondition = Utils.dget('byCondition', self.properties, self.byCondition)
        self.removeOriginal = Utils.dget('remove', self.properties, self.removeOriginal)
        return True

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nmerge = 0
        for smp in SC.samples:
            smp['bam'] = smp['name'] + ".bam"
            LOG.log("Concatenating BAM files for sample {} into {}", smp['name'], smp['bam'])
            cmdline = "generic.qsub samtools cat -o " + smp['bam']
            for rs in smp['readsets']:
                cmdline += " " + rs['bam']
            cmdline += " module:samtools"
            if not self.dry:
                ACT.submit(cmdline, done='bamcat.@.done')
                nmerge += 1
        ACT.wait(('bamcat.@.done', nmerge))

        if self.byCondition:
            nmerge = 0
            for cond in SC.conditions:
                cond['bam'] = cond['name'] + ".bam"
                LOG.log("Concatenating BAM files for condition {} into {}", cond['name'], cond['bam'])
                cmdline = "generic.qsub samtools cat -o " + cond['bam']
                for sb in SC.conditionBAMs(cond['name']):
                    cmdline += " " + sb
                cmdline += " module:bamtools"
                if not self.dry:
                    ACT.submit(cmdline, done='cbamcat.@.done')
                    nmerge += 1
            ACT.wait(('cbamcat.@.done', nmerge))

        return True

### RSEM 

class RSEMquant(Line):
    """Perform transcriptome quantification using RSEM. Calls rsem.qsub on the BAM files for each sample."""
    name = "RSEM (quantification)"
    version = ""

    def Verify(self):
        ACT = self.actor
        version = ACT.shell("module load rsem; rsem-calculate-expression --version")
        if version.startswith("Current version:"):
            self.version = version[17:]
            return True
        else:
            return False

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nrsem = 0
        for smp in SC.samples:
            smp['rsembam'] = smp['name'] + ".transcript.bam"
            paired = smp['readsets'][0]['paired'] # we assume that all readsets in a sample are the same
            if paired:
                cmdline = "rsem.qsub {} {} bamp={}".format(smp['name'], ACT.rsemidx, smp['bam'])
            else:
                cmdline = "rsem.qsub {} {} bam={}".format(smp['name'], ACT.rsemidx, smp['bam'])
            if not self.dry:
                ACT.submit(cmdline, done="rsem.@.done")
                nrsem += 1
        return ACT.wait(("rsem.@.done", nrsem))

class RSEMdiff(Line):
    """Perform differential analysis using RSEM. """
    name = "RSEM (differential)"
    version = ""

    def Verify(self):
        ACT = self.actor
        version = ACT.shell("module load rsem; rsem-calculate-expression --version")
        if version.startswith("Current version:"):
            self.version = version[17:]
            return True
        else:
            return False

    def writeRSEMqsubs(self, ACT, contrast):
        SC = ACT.sc
        ctrlcond = contrast['control']
        ctrlSamples = SC.conditionSamples(ctrlcond)
        testcond = contrast['test']
        testSamples = SC.conditionSamples(testcond)
        nctrl = len(ctrlSamples)
        ntest = len(testSamples)
        contrast['rsemgqsub'] = gfilename = "{}.vs.{}.g.qsub".format(testcond, ctrlcond)
        contrast['gmatrix'] = "{}.vs.{}.gmatrix.csv".format(testcond, ctrlcond)
        contrast['gdiff'] = "{}.vs.{}.gdiff.csv".format(testcond, ctrlcond)
        contrast['gfdr'] = "{}.vs.{}.gfdr.csv".format(testcond, ctrlcond)

        with open(gfilename, "w") as out:
            out.write("""#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=10G

module load rsem/1.2.29
module load dibig_tools

#rsem-generate-data-matrix 

rnasequtils.py matrix """)
            if ACT.erccdb:
                out.write(" -e -ercc " + ACT.erccdb + " -mix " + ",".join([s['mix'] for s in ctrlSamples + testSamples]))
            for s in ctrlSamples + testSamples:
                out.write(" {}.genes.results".format(s['name']))
            out.write(" > {}\n\n".format(contrast['gmatrix']))
            out.write("rsem-run-ebseq {} {},{} {}\n".format(contrast['gmatrix'], nctrl, ntest, contrast['gdiff']))
            out.write("rsem-control-fdr {} {} {}\n".format(contrast['gdiff'], ACT.fdr, contrast['gfdr']))

        contrast['rsemiqsub'] = ifilename = "{}.vs.{}.i.qsub".format(testcond, ctrlcond)
        contrast['imatrix'] = "{}.vs.{}.imatrix.csv".format(testcond, ctrlcond)
        contrast['idiff'] = "{}.vs.{}.idiff.csv".format(testcond, ctrlcond)
        contrast['ifdr'] = "{}.vs.{}.ifdr.csv".format(testcond, ctrlcond)

        with open(ifilename, "w") as out:
            out.write("""#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=10G

module load rsem/1.2.29
module load dibig_tools

#rsem-generate-data-matrix 

rnasequtils.py matrix """)
            if ACT.erccdb:
                out.write(" -e -ercc " + ACT.erccdb + " -mix " + ",".join([s['mix'] for s in ctrlSamples + testSamples]))
            for s in ctrlSamples + testSamples:
                out.write(" {}.isoforms.results".format(s['name']))
            out.write(" > {}\n\n".format(contrast['imatrix']))
            out.write("rsem-run-ebseq {} {},{} {}\n".format(contrast['imatrix'], nctrl, ntest, contrast['idiff']))
            out.write("rsem-control-fdr {} {} {}\n".format(contrast['idiff'], ACT.fdr, contrast['ifdr']))
        return (gfilename, ifilename)

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nrsem = 0
        for contr in SC.contrasts:
            (gqsub, iqsub) = self.writeRSEMqsubs(ACT, contr)
            if not self.dry:
                ACT.submit(gqsub, done='rsem.@.done')
                ACT.submit(iqsub, done='rsem.@.done')
                nrsem += 2
        return ACT.wait(("rsem.@.done", nrsem))

    def filterDiff(self, infile, outfile, fc, translation=None, extra=[], biotype=None):
        """Read RSEM differential expression results from `infile' and write them to `outfile', converting
the fold changes to log2 and discarding records having log2(FC) < `fc'. If `translation' is provided
it should be a dictionary mapping gene or transcript IDs (first column) to additional annotations. The
`extra' list provides the header names for the additional annotations. If `biotype' is specified, the
annotations are assumed to contain a 'biotype' as the second element, and the output file will contain only records having
the specified biotype. Returns a tuple containing the number of records written, number of overexpressed,
number of underexpressed, and a dictionary mapping each biotype to is number of occurrences."""
        data = []
        biotypes = {}
        missing = ["???" for r in range(len(extra)) ]
        ng = nup = ndown = 0
        with open(infile, "r") as f:
            f.readline()
            for line in f:
                good = True
                parsed = line.rstrip("\r\n").split("\t")
                log2fc = math.log(float(parsed[3]), 2)
                if abs(log2fc) >= fc:
                    if translation:
                        k = parsed[0].strip('"')
                        if k in translation:
                            trans = translation[k]
                            biot = trans[1]
                            if biotype:
                                if biotype != biot:
                                    good = False
                                else:
                                    Utils.dinc(biot, biotypes)
                        else:
                            trans = missing
                        if good:
                            row = [log2fc] + trans + [str(log2fc), parsed[1], parsed[5], parsed[6]]
                    else:
                        row = [log2fc, parsed[0], str(log2fc), parsed[1], parsed[5], parsed[6]]

                    if good:
                        ng += 1
                        if log2fc > 0:
                            nup += 1
                        else:
                            ndown += 1
                        data.append(row)

        data.sort(key=lambda d: d[0], reverse=True)

        with open(outfile, "w") as out:
            if extra:
                out.write("\t".join(extra) + "\tLog2(FC)\tP-value\tExp0\tExp1\n")
            else:
                out.write("Gene\tLog2(FC)\tP-value\tExp0\tExp1\n")
            for d in data:
                out.write("\t".join(d[1:]) + "\n")
        return (ng, nup, ndown, biotypes)

    def PostExecute(self):
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log

        LOG.log("Parsing GTF file {}", ACT.gtf)
        (gtable, txtable) = Utils.parseGTF(ACT.gtf)
        allfinalg = []
        allfinalc = []
        allfinali = []
        mergecmdline = "module load dibig_tools; rnasequtils.py merge "
        LOG.log("Generating final differential expression tables.")
        for contr in SC.contrasts:
            ctrlcond = contr['control']
            testcond = contr['test']
            vs = "{}.vs.{}".format(testcond, ctrlcond)
            contr['genefinal'] =   vs + ".geneDiff.csv"
            contr['codingfinal'] = vs + ".codinggeneDiff.csv"
            contr['isofinal'] =    vs + ".isoDiff.csv"
            mergecmdline += vs + " "
            allfinalg.append(contr['genefinal'])
            allfinalc.append(contr['codingfinal'])
            allfinali.append(contr['isofinal'])
            contr['ngdiff'] = self.filterDiff(contr['gfdr'], contr['genefinal'], ACT.fc, translation=gtable, extra=["ENSG", "Biotype", "Gene"])
            contr['ncdiff'] = self.filterDiff(contr['gfdr'], contr['codingfinal'], ACT.fc, translation=gtable, extra=["ENSG", "Biotype", "Gene"],
                                              biotype='protein_coding')
            contr['nidiff'] = self.filterDiff(contr['ifdr'], contr['isofinal'], ACT.fc, translation=txtable, extra=["ENST", "Biotype", "Transcript", "Gene"])
        ACT.shell(mergecmdline)
        ACT.shell("module load dibig_tools; csvtoxls.py merged.allDiff.xlsx -q merged.geneDiff.csv merged.codinggeneDiff.csv merged.isoDiff.csv")
        ACT.shell("module load dibig_tools; csvtoxls.py genediff.xlsx -q {}".format(" ".join(allfinalg)))
        ACT.shell("module load dibig_tools; csvtoxls.py codingdiff.xlsx -q {}".format(" ".join(allfinalc)))
        ACT.shell("module load dibig_tools; csvtoxls.py isodiff.xlsx -q {}".format(" ".join(allfinali)))
        return True

    def Report(self):
        ACT = self.actor
        SC = ACT.sc

        ACT.scene("Differential expression - gene level")
        ACT.reportf("""Differential gene expression was analyzed using <b>{}</b>. The following table reports the number of differentially expressed 
genes in each contrast with <b>abs(log<inf>2</inf>(FC)) >= {}</b> and <b>FDR-corrected P-value <= {}</b>.""".format(self.version, ACT.fc, ACT.fdr))
        ACT.reportf("""The lists of differentially expressed genes for all contrasts can also be downloaded as a single Excel file using the link below.""")
        tbl4 = Table.ScrollingTable(id='tbl4', align="LLRRRCC", caption="Results of gene-level differential expression analysis.")
        tbl4.startHead()
        tbl4.addHeaderRow(["Test", "Control", "Overexpressed", "Underexpressed", "Total", "Table", "Expressions"])
        tbl4.startBody()
        for contr in SC.contrasts:
            ng = contr['ncdiff']
            tbl4.addRow([ contr['test'], 
                          contr['control'], 
                          Utils.UP(ng[1]), Utils.DOWN(ng[2]), ng[0],
                          Utils.linkify(contr['codingfinal']),
                          Utils.linkify(contr['gmatrix'])])
        tbl4.toHTML(ACT.out)
        ACT.file("codingdiff.xlsx", description="Excel file containing differentially expressed genes for all contrasts (one sheet per contrast). Only includes coding genes.")

        ACT.reportf("""The following table reports results from the same differential analysis as above, but includes all biotypes instead of coding genes only.""")

        tbl4b = Table.ScrollingTable(id='tbl4b', align="LLRRRCC", caption="Results of gene-level differential expression analysis (all biotypes).")
        tbl4b.startHead()
        tbl4b.addHeaderRow(["Test", "Control", "Overexpressed", "Underexpressed", "Total", "Table", "Expressions"])
        tbl4b.startBody()
        for contr in SC.contrasts:
            ng = contr['ngdiff']
            tbl4b.addRow([ contr['test'], 
                           contr['control'], 
                           Utils.UP(ng[1]), Utils.DOWN(ng[2]), ng[0],
                           Utils.linkify(contr['genefinal']),
                           Utils.linkify(contr['gmatrix'])])
        tbl4b.toHTML(ACT.out)
        ACT.file("genediff.xlsx", description="Excel file containing differentially expressed genes for all contrasts (one sheet per contrast). Includes all genes and pseudo-genes.")

        ACT.scene("Differential expression - isoform level")
        ACT.reportf("""Differential isoform expression was analyzed using <b>{}</b>. The following table reports the number of differentially expressed 
isoforms in each contrast with <b>abs(log<inf>2</inf>(FC)) >= {}</b> and <b>FDR-corrected P-value <= {}</b>.""".format(self.version, ACT.fc, ACT.fdr))
        ACT.reportf("""The lists of differentially expressed isoforms for all contrasts can also be downloaded as a single Excel file using the link below.""")
        tbl5 = Table.ScrollingTable(id='tbl5', align="LLRRRCC", caption="Results of isoform-level differential expression analysis.")
        tbl5.startHead()
        tbl5.addHeaderRow(["Test", "Control", "Overexpressed", "Underexpressed", "Tot isoforms", "Table", "Expressions"])
        tbl5.startBody()
        for contr in SC.contrasts:
            ni = contr['nidiff']
            tbl5.addRow([ contr['test'], 
                          contr['control'], 
                          Utils.UP(ni[1]), Utils.DOWN(ni[2]), ni[0],
                          Utils.linkify(contr['isofinal']),
                          Utils.linkify(contr['imatrix']) ])
        tbl5.toHTML(ACT.out)
        ACT.file("isodiff.xlsx", description="Excel file containing differentially expressed isoforms for all contrasts (one sheet per contrast).")
        return True

### BAMtoWig

class BAMtoWig(Line):
    """Convert BAM files for samples to Wig files."""
    name = "BAM to Wig"
    scale = 1000000000
    countsfile = "afterCounts.csv"

    def Setup(self):
        self.scale = int(Utils.dget('scale', self.properties, self.scale))
        self.countsfile = Utils.dget('countsfile', self.properties, self.countsfile)
        return True

    def Execute(self):
        return True

    def PostExecute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        counts = Utils.fileToDict(self.countsfile)
        nwig = 0
        for smp in SC.samples:
            name = smp['name']
            smp['wigfile'] = name + ".wig"
            if not self.dry:
                ACT.submit("""generic.qsub bamToWig.py -o {} -n {} -s {} -t {} -d "{} transcriptome" {} module:dibig_tools""".format(
                        smp['wigfile'], counts[name], self.scale, name, name, smp['genomebam']), done="wig.@.done")
                nwig += 1
        return ACT.wait(("wig.@.done", nwig))

### MACS

class MACScallpeak(Line):
    """Run MACS in callpeak mode."""
    name = "MACS (callpeak)"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        # print ACT.sc.contrasts
        # print ACT.sc.conditions
        # print ACT.sc.samples
        # raw_input()

        nmacs = 0
        for cond in SC.conditions:
            cond['macsdir'] = macsdir = cond["name"] + ".macs"
            bams = ",".join(SC.conditionBAMs(cond['name']))
            if not self.dry:
                ACT.submit("macs.qsub {} {} {}".format(cond['name'], macsdir, bams), done="macs.@.done")
                nmacs += 1

        for contrast in SC.contrasts:
            ctrlcond = contrast['control']
            testcond = contrast['test']
            name = "{}.vs.{}".format(testcond, ctrlcond)
            contrast['macsdir'] = macsdir = "{}.vs.{}.macs".format(testcond, ctrlcond)
            cbams = ",".join(SC.conditionBAMs(ctrlcond))
            tbams = ",".join(SC.conditionBAMs(testcond))
            if not self.dry:
                ACT.submit("macs-diff.qsub {} {} {} {}".format(name, macsdir, tbams, cbams), done="macs.@.done")
                nmacs += 1
        return ACT.wait(("macs.@.done", nmacs))
    
### Methylation

class CSfilter(Line):
    """Invoke cscall -filter on all fastq files. This call replaces the 'left' and 'right'
fastq files with the filtered fastqs, while the original ones are stored under 'leftprefilter'
and 'rightprefilter'."""
    name = "CSCALL filter"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nfilter = 0
        if ACT.sites:
            sites = " ".join(ACT.sites.split(","))
            for rs in ACT.sc.readsets:
                if rs['paired']:
                    # Paired-end
                    rs['leftprefilter'] = rs['left']
                    rs['rightprefilter'] = rs['right']
                    rs['left'] = "good-" + os.path.split(rs['left'])[1]
                    rs['right'] = "good-" + os.path.split(rs['right'])[1]
                    LOG.log("Filtering unconverted reads from `{}' and `{}' for sites {}".format(rs['leftprefilter'], rs['rightprefilter'], sites))
                    if not self.dry:
                        ACT.submit("generic.qsub cscall -filter -1 {} -2 {} -amax {} -fmax {} -s {} module:dibig_tools".format(
                                rs['leftprefilter'], rs['rightprefilter'], ACT.amax, ACT.fmax, sites), done="filter.@.done")
                        nfilter += 1
                else:
                    # Single-end
                    rs['leftprefilter'] = rs['left']
                    rs['left'] = "good-" + os.path.split(rs['left'])[1]
                    LOG.log("Filtering unconverted reads from `{}' for sites {}".format(rs['leftprefilter'], sites))
                    if not self.dry:
                        ACT.submit("generic.qsub cscall -filter -1 {} -amax {} -fmax {} -s {} module:dibig_tools".format(
                                rs['leftprefilter'], ACT.amax, ACT.fmax, sites), done="filter.@.done")
                        nfilter += 1
            return ACT.wait(("filter.@.done", nfilter))
        else:
            return True

    def Report(self):
        ACT = self.actor
        SC = ACT.sc

        ACT.scene("Read pre-filtering")
        ACT.reportf("""The input sequences were filtered using <b>CSCALL</b> to remove reads containing more than <b>{}</b> or more than <b>{}%</b> 'lone' cytosines.
Lone cytosines are defined as cytosines that are not part in one of the following sites: <b>{}</b>.""".format(ACT.amax, ACT.fmax*100, ACT.sites))

        prereadcounts = Utils.fileToDict(ACT.fastqCounts)
        readcounts = Utils.fileToDict(ACT.fastqCountsFilter)
        tbl2 = Table.ScrollingTable(id='tblcsfilter', align="LRRRR", caption='Number of reads before and after cscall filtering.')
        tbl2.startHead()
        tbl2.addHeaderRow(["Sample", "Reads before filtering", "Reads after filtering", "Reads removed", "Pct retained"])
        tbl2.startBody()
        for smp in SC.samples:
            sn = smp['name']
            tbl2.addRow([ "<b>" + sn + "</b>", 
                          Utils.fmt(prereadcounts[sn]),
                          Utils.fmt(readcounts[sn]), 
                          Utils.fmt(prereadcounts[sn] - readcounts[sn]),
                          Utils.pct(readcounts[sn], prereadcounts[sn])])
        tbl2.toHTML(ACT.out)
        return True

### MMAP

class MMAP(Line):
    """Run BSMAP through MMAP."""
    name = "BSMAP"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc
        
        nmap = 0
        writeBSMAPscripts(ACT)

        for smp in SC.samples:
            smp['bam'] = smp['name'] + ".bam"
            LOG.log("Calling BSMAP on `{}'.", smp['bam'])
            if not self.dry:
                j=ACT.submit(smp['bsmapscript'])
                ACT.submit("bam-sort-and-index.qsub {} Y".format(smp['bam']), after=j, done="mmap.@.done")
                nmap += 1
        return ACT.wait(("mmap.@.done", nmap))

    def Report(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nalign = Utils.fileToDict(ACT.alignedCounts)
        readcounts = Utils.fileToDict(ACT.fastqCountsFilter)

        ACT.scene("Short-read alignment")
        ACT.reportf("""The input reads were aligned to the genome using <b>bsmap</b>. The following table reports the number of total and aligned reads for each sample.""")
        tbl = Table.ScrollingTable(id='tblbsmap', align="LRRR", caption="Number of input and aligned reads, with percentage of aligned reads.")
        tbl.startHead()
        tbl.addHeaderRow(["Sample", "Total reads", "Aligned reads", "Pct Aligned"])
        tbl.startBody()
        for smp in SC.samples:
            sn = smp['name']
            tbl.addRow([ "<b>" + sn + "</b>",
                         Utils.fmt(readcounts[sn]*2),                # *** HACK! assuming paired-end
                         Utils.fmt(nalign[sn]),
                         Utils.pct(nalign[sn], readcounts[sn]*2) ])
        tbl.toHTML(ACT.out)
        return True
        
### CScall

class CScall(Line):
    """Perform methylation calling with CSCALL."""
    name = "CSCALL"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        ncall = 0
        for cond in SC.conditions:
            bams = ",".join([ smp['bam'] for smp in SC.conditionSamples(cond['name']) ])
            cond['callbed'] = cond['name'] + ".bed"
            cond['callreport'] = cond['name'] + "-report.csv"
            cond['lcreport'] = cond['name'] + "-lc.csv"
            cond['histreport'] = cond['name'] + "-hist.csv"
            if ACT.strand == None:
                strand = ""
            else:
                strand = "strand=" + ACT.strand
            LOG.log("Calling methylation on BAM files for sample {}.".format(cond['name']))
            if not self.dry:
                ACT.submit("cscall.qsub {} {} {} out={} report={} lc={} hist={} mind={} mins={} {}".format(bams, ACT.siteindex, ACT.reference, cond['callbed'], cond['callreport'], cond['lcreport'], cond['histreport'], ACT.mindepth, ACT.minsamples, strand), done="mcall.@.done")
                ncall += 1
        return ACT.wait(("mcall.@.done", ncall))

    def PostExecute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        for cond in SC.conditions:
            cond['callreportx'] = cond['name'] + "-report.xlsx"
            ACT.shell("module load dibig_tools; csvtoxls.py {} {} -firstrowhdr -name Report {} -firstrowhdr -name Histogram {} -firstrowhdr -name LoneC",
                      cond['callreportx'], cond['callreport'], cond['histreport'], cond['lcreport'])

        return True

    def Report(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc
        
        ACT.scene("Methylation calling")
        ACT.reportf("""Methylation calling was performed with <b>cscall v1.0</b> using the <b>{}</b> index. 
A site was included in the analysis if it reached a coverage of <b>{}</b> in at least <b>{}</b> sample{}.""", 
                    os.path.splitext(os.path.basename(ACT.siteindex))[0],
                    ACT.mindepth, ACT.minsamples, "" if ACT.minsamples == 1 else "s")
        ACT.reportf("""The following table provides a summary of the number of sites identified for each experimental condition. Follow the links to the reports for more detailed information.""")
        tbx = Table.ScrollingTable(id='tblcscall', align="LRRRRCCCC", caption="Summary of methylation calling.")
        tbx.startHead()
        tbx.addHeaderRow(["Condition", "Sites", "Total basepairs", "Total site coverage", "Avg site coverage", "Full report", "Plots"])
        tbx.startBody()
        for cond in SC.conditions:
            data = getCscallReport(cond)
            if data:
                tbx.addRow([ "<b>" + cond['name'] + "</b>",
                             Utils.fmt(int(data[1])),
                             Utils.fmt(int(data[0])),
                             Utils.fmt(int(data[2])),
                             data[3],
                             Utils.linkify(cond['callreportx']),
                             Utils.linkify(cond['methplots']) ])
        tbx.toHTML(ACT.out)

        # ACT.reportf("""The <i>Sites</i> column contains the number of sites for the desired methylation type found in the alignment. """)
        return True
        
### MCOMP

class MCOMP(Line):
    name = "MCOMP"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        ncomp = 0
        subs = writeMCOMPqsub(ACT)
        if not self.dry:
            for sub in subs:
                ncomp += 1
                ACT.submit(sub, done="mcomp.@.done")
        ACT.wait(("mcomp.@.done", ncomp))
        LOG.log("Computing differential methylation summary.")
        self.prefix = Utils.id_generator(10)
        ACT.diffmeth = diffmethSummary(ACT, self.prefix, self.dry)
        return True

    def Report(self):
        ACT = self.actor
        LOG = ACT.log

        ACT.scene("Differential methylation analysis")
        ACT.reportf("""Differential methylation was determined using the <b>mcomp</b> program. Sites for which the P-value of the difference between test and control methylation rates was below <b>{}</b> were considered significant. Sites with a methylation difference above <b>{}</b> were classified as <i>Strongly hyper/hypo methylated</i> (++ or --), otherwise they were classified as <i>hyper/hypo methylated</i> (+ or -). The following table reports the number of differentially methylated sites identified in each contrast, and provides links to full reports.""".format(ACT.pval, ACT.diff))
        tbl = Table.ScrollingTable(id='tblmcomp', align="CCRRRRRR", caption="Results of differential methylation analysis. The <i>Total sites</i> column contains the total number of sites examined for each sample, while <i>Significant sites</i> shows how many sites display evidence of statistically significant differential methylation. ++ = strong hypermethylation, + = hypermethylation, - = hypomethylation, -- = strong hypomethylation.")
        tbl.startHead()
        tbl.addHeaderRow(["Test", "Control", "Total Sites", "Significant sites", "++", "+", "-", "--"])
        tbl.startBody()
        for d in ACT.diffmeth:
            totsig = d['++'] + d['+'] + d['-'] + d['--']
            tbl.addRow([ "<b>" + d['test'] + "</b>",
                         "<b>" + d['control'] + "</b>",
                         Utils.fmt(d['total']),
                         "{} ({})".format(totsig, Utils.pct(totsig, d['total'])),
                         "{} ({})".format(Utils.UP(Utils.fmt(d['++'])), Utils.pct(d['++'], totsig)),
                         "{} ({})".format(Utils.up(Utils.fmt(d['+'])), Utils.pct(d['+'], totsig)),
                         "{} ({})".format(Utils.down(Utils.fmt(d['-'])), Utils.pct(d['-'], totsig)),
                         "{} ({})".format(Utils.DOWN(Utils.fmt(d['--'])), Utils.pct(d['--'], totsig))
                         ])

        tbl.toHTML(ACT.out)

        ACT.reportf("""<br><br>The following table provides links to files containing detailed differential methylation data. The <i>Significant</i> files contain sites showing differential methylation only, while the <i>Full</i> files contain data for all sites. The <i>BED</i> files can be loaded into a genome browser (e.g., IGV), and the <i>Plots</i> files contain plots of differential methylation across chromosomes.""")

        tbl = Table.ScrollingTable(id='tblmcomp2', align="CCCCCC", caption="Tables of differential methylation data.")
        tbl.startHead()
        tbl.addHeaderRow(["Test", "Control", "Significant", "Full", "BED", "Plots"])
        tbl.startBody()
        for d in ACT.diffmeth:
            tbl.addRow([ "<b>" + d['test'] + "</b>",
                         "<b>" + d['control'] + "</b>",
                         Utils.linkify(d['sigfilex']),
                         Utils.linkify(d['fullfilex']),
                         Utils.linkify(d['bedfile']),
                         Utils.linkify(d['dmcfile']) ])

        tbl.toHTML(ACT.out)

        return True

### GeneMeth

class GeneMeth(Line):
    name = "Gene methylation"

    def makeCmdline(self, ACT, infile, outfile, srt=False):
        cmdline = "generic.qsub genediffmeth "
        if ACT.regions:
            cmdline += " -r " + ACT.regions
        if ACT.size:
            cmdline += " -d " + ACT.size
        if ACT.mode:
            cmdline += " -m " + ACT.mode
        if ACT.nsites:
            cmdline += " -l " + ACT.nsites
        if srt:
            cmdline += " -s "
        return cmdline + " " + infile + " " + ACT.genesdb + " " + outfile + " module:dibig_tools"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        
        ngd = 0

        if ACT.genesdb:
            for dm in ACT.diffmeth:
                outfile = "{}.vs.{}.genes.sig.csv".format(dm['test'], dm['control'])
                dm['genesigfile'] = outfile
                cmdline = self.makeCmdline(ACT, dm['sigfile'], outfile, True)
                LOG.log("Generating differential methylation file for genes containing significant sites.")
                if not self.dry:
                    ACT.submit(cmdline, done="gdm.@.done")
                    ngd += 1
                outfile = "{}.vs.{}.genes.full.csv".format(dm['test'], dm['control'])
                dm['genefullfile'] = outfile
                cmdline = self.makeCmdline(ACT, dm['fullfile'], outfile)
                LOG.log("Generating differential methylation file for all genes.")
                if not self.dry:
                    ACT.submit(cmdline, done="gdm.@.done")
                    ngd += 1
        return ACT.wait(("gdm.@.done", ngd))

    def PostExecute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        for dm in ACT.diffmeth:
            dm['genesigfilex'] = "{}.vs.{}.genes.sig.xlsx".format(dm['test'], dm['control'])
            dm['genefullfilex'] = "{}.vs.{}.genes.full.xlsx".format(dm['test'], dm['control'])
            ACT.shell("module load dibig_tools; csvtoxls.py {} -q {}; csvtoxls.py {} -q {}",
                      dm['genesigfilex'], dm['genesigfile'], 
                      dm['genefullfilex'], dm['genefullfile'])
        return True
        
    def countSigGenes(self, filename, diff):
        counts = {'++': 0, '+': 0, '-': 0, '--': 0, 'total': 0}
        with open(filename, "r") as f:
            f.readline()
            for line in f:
                line = line.rstrip("\n\r").split("\t")
                d = float(line[6])
                if d > 0:
                    if d > diff:
                        tag = '++'
                    else:
                        tag = '+'
                else:
                    if d < -diff:
                        tag = '--'
                    else:
                        tag = '-'
                counts['total'] += 1
                counts[tag] += 1
        return counts

    def Report(self):
        ACT = self.actor
        LOG = ACT.log

        if ACT.genesdb:
            ACT.scene("Gene methylation analysis")
            ACT.reportf("""The following table provides links to files containing information about methylation in genes. For each contrast, two files are provided, one 
    containing only genes showing evidence of differential methylation (ordered from highest hypermethylation to highest hypomethylation) and one containing all
    genes (ordered by chromosome and position). Gene differential methylation was computed collecting all differentially methylated sites in the <b>{}</b> region
of the gene (where p = promoter, b = gene body, d = downstream) and combining them with the <b>{}</b> operator (only including genes containing at least <b>{}</b> 
site{}). Promoter and downstream were defined as <b>{}bp</b> upstream and downstream of the transcript, respectively.""".format(
                    ACT.regions, ACT.mode, ACT.nsites, Utils.plural(ACT.nsites), ACT.size))
            tbl = Table.ScrollingTable(id='tblgmeth', align="CCRRRRRCC", caption="Differential methylation in genes. For each contrast, the table lists the number of genes showing evidence of differential methylation, and the breakdown into the four differential methylation classes.")
            tbl.startHead()
            tbl.addHeaderRow(["Test", "Control", "Genes", "++", "+", "-", "--", "Significant", "Full"])
            tbl.startBody()
            for d in ACT.diffmeth:
                counts = self.countSigGenes(d['genesigfile'], ACT.diff)
                print counts
                tot = counts['total']
                tbl.addRow([ "<b>" + d['test'] + "</b>",
                             "<b>" + d['control'] + "</b>",
                             Utils.fmt(tot),
                             "{} ({})".format(Utils.UP(Utils.fmt(counts['++'])),   Utils.pct(counts['++'], tot)),
                             "{} ({})".format(Utils.up(Utils.fmt(counts['+'])),    Utils.pct(counts['+'], tot)),
                             "{} ({})".format(Utils.down(Utils.fmt(counts['-'])),  Utils.pct(counts['-'], tot)),
                             "{} ({})".format(Utils.DOWN(Utils.fmt(counts['--'])), Utils.pct(counts['--'], tot)),
                             Utils.linkify(d['genesigfilex']),
                             Utils.linkify(d['genefullfilex']) ])

            tbl.toHTML(ACT.out)
        return True

### Methylation plots

class MethPlots(Line):
    name = "Plots of absolute methylation level."

    def Execute(self):
        ACT = self.actor
        SC = ACT.sc

        for cond in SC.conditions:
            drawBEDplots(ACT, cond, self.dry)

        drawDMCplots(ACT, self.dry)
        return True

### For BA3P

class ContigRenamer(Line):
    """Rename sequences in contig files as XXX_n, where XXX is the sample name and n is a progressive number."""
    name = "Rename contigs"

    def Execute(self):
        ACT = self.actor
        SC = ACT.sc

        ACT.shell("touch dummymap")
        for smp in SC.samples:
            infasta = smp['fasta']
            outfasta = "Contigs/" + smp['name'] + ".spades.fasta"
            print outfasta
            smp['contig'] = outfasta
            if not self.dry:
                ACT.shell("module load dibig_tools; fastools -map -mapfile-in dummymap -in {} -out {} -if-missing counter -prefix {}_", infasta, outfasta, smp['name'])
        return True

class IndexReference(Line):
    """This line checks if the reference file in ACT.reference needs to be indexed with samtools."""
    name = "Index reference"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        fai = ACT.reference + ".fai" 
        refdict = ACT.setFileExt(ACT.reference, ".dict")

        if ACT.missingOrStale(fai, ACT.reference):
            LOG.log("Reference index {} does not exist or is out of date, creating it.", fai)
            if not self.dry:
                ACT.shell("module load samtools; samtools faidx {}", ACT.reference)

        if ACT.missingOrStale(refdict, ACT.reference):
            LOG.log("Reference dictionary {} does not exist, creating it.", refdict)
            if not self.dry:
                ACT.submit("picard.qsub CreateSequenceDictionary R={} O={}".format(ACT.reference, refdict), done="dict.done")
                ACT.wait("dict.done")
        return True

class Prokka(Line):
    name = "Prokka"

    def Verify(self):
        ACT = self.actor

        if ACT.genus == None or ACT.species == None:
            self.status = "prokka requires specification of genus and species in configuration file."
            return False
        else:
            return True
        
    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        ACT.nprokka = 0
        for smp in SC.samples:
            prokkadir = smp['name']
            smp['prokkadir'] = prokkadir
            fasta = smp['name'] + ".spades.fasta"
            if ACT.locus == None:
                locus = ""
            else:
                locus = "locus=" + ACT.locus
            if not self.dry:
                os.chdir("Contigs")
                ACT.submit('prokka.qsub {} {} "{}" "{}" {} prefix={}'.format(fasta, prokkadir, ACT.genus, ACT.species, locus, smp['name']), done="../prokka.@.done")
                os.chdir("..")
                ACT.nprokka += 1
        return True

class Mauve(Line):
    """Run mauve contig mover on all fasta files in the Contigs directory."""
    name = "Mauve contig mover"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        LOG.log("Running Mauve contig mover on directory Contigs/")
        nmauve = 0
        if not self.dry:
            # ACT.submit("mauve-contig-mover.qsub {}".format(ACT.reference), done="mauve.done")
            # ACT.wait("mauve.done")
            for smp in SC.samples:
                ACT.submit("mauve-contig-mover.qsub {} {}".format(ACT.reference, smp['contig']), done="mauve.@.done")
                nmauve += 1
            ACT.wait(("mauve.@.done", nmauve))
        LOG.log("Selecting optimal alignment for each sample.")
        ACT.mkdir("Alignments")
        for smp in SC.samples:
            base = "Ordered/" + smp['name'] + ".spades/"
            LOG.log("{}: examining directory {}.", smp['name'], base)
            alns = os.listdir(base)
            maxa = [0, '']
            for a in alns:
                if a[0:9] == 'alignment':
                    i = int(a[9:])
                    if i > maxa[0]:
                        maxa = [i, a]
            LOG.log("{}: best alignment is {}", smp['name'], maxa[0])
            ACT.shell("cp {}/{}/{}.spades.fasta Alignments/".format(base, maxa[1], smp['name']))

        return True
    
class ProgMauve(Line):
    """Run progressive-mauve on all files in the Alignments directory. Files listed in the 'additional'
property are added to the directory before running progressive mauve."""
    name = "Progressive Mauve"
    
    def Execute(self):
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log
        
        if not self.dry:
            os.chdir("Alignments")
            if ACT.additional:
                for a in ACT.additional:
                    if os.path.exists(a):
                        LOG.log("Copying additional reference {} into Alignments.", a)
                        ACT.shell("cp {} .".format(a))
                    else:
                        LOG.log("Additional reference {} not found.", a)
            ACT.submit("progressive-mauve.qsub denovo-final denovo-final", done="../progr.done")
            os.chdir("..")
            ACT.wait("progr.done")
        return True

    def PostExecute(self):
        ACT = self.actor
        
        ## Convert final xmfa to fasta
        ACT.shell("module load dibig_tools perl; convertAlignment.pl -i Alignments/denovo-final.xmfa -o denovo-final.fasta -f fasta -g xmfa -c")
        return True

class Roary(Line):
    name = "Roary"
    
    def Execute(self):
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log

        ACT.wait(("prokka.@.done", ACT.nprokka))
        if ACT.nprokka > 0:
            ACT.mkdir("Roary")
            for smp in SC.samples:
                ACT.shell("cp {}/*.gff Roary/".format(smp['prokkadir']))
            if not self.dry:
                os.chdir("Roary")
                ACT.submit("roary.qsub", done="../roary.done")
                os.chdir("..")
                ACT.wait("roary.done")

            ACT.shell("cp Roary/core_gene_alignment.aln .")
            ACT.shell("module load snp-sites; snp-sites -m -o core-alignment-snps.fa core_gene_alignment.aln")
        return True

class SNPsites(Line):
    name = "snp-sites"
    infile = None
    outfile = None

    def Setup(self):
        self.infile = Utils.dget('infile', self.properties, self.infile)
        self.outfile = Utils.dget('outfile', self.properties, self.outfile)

    def Verify(self):
        if self.infile and self.outfile:
            return True
        self.status = "Need to specify infile and outfile!"
        return False

    def Execute(self):
        ACT = self.actor

        ACT.shell("module load snp-sites; snp-sites -m -o {} {}".format(self.outfile, self.infile))
        return True
