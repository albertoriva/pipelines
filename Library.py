# (c) 2015, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import os
import math

import Utils
import PipelineUtils
from Lines import Line
from SampleCollection import SampleCollection
from dmaputils import initializeContrasts, writeBSMAPscripts, writeMCOMPqsub, diffmethSummary, drawBEDplots, drawDMCplots, getCscallReport, runAvgMeth, runMethHist
from rnasequtils import writeMatrixScript, writeContrastMatrix, filterDiff, exprScatterplot
from chipsequtils import readHomerAnnots, writeClassifyScript, classificationBarChart
from snpcallutils import splitVCF, snpeffRec
import Table
import GeneList
import Intersector

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
        LOG = ACT.log
        SC = ACT.sc

        for smp in SC.samples:
            smp['mix'] = ACT.getConf('mix', section=smp['name'], default='1')
        LOG.log("Loading genes from database {}", ACT.genesdb)
        loader = GeneList.getLoader(ACT.genesdb)
        if loader:
            ld = loader(ACT.genesdb)
            ACT.genelist = ld.load()
            LOG.log("{} genes loaded.", ACT.genelist.ngenes)
            ACT.gtable  = ACT.genelist.getGenesTable()
            ACT.txtable = ACT.genelist.getTranscriptsTable()
            return True
        else:
            return False
    
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
        self.adapter = Utils.dget('adapter', self.properties, None)
        self.minlen = Utils.dget('minlen', self.properties, None)
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
        if self.adapter:
            adapt = "ad=" + self.adapter
        else:
            adapt = ""
        if self.minlen:
            minlen = "minlen=" + self.minlen
        else:
            minlen = ""

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
                    rs['fastqcbeforeleft'] = PipelineUtils.fastqcPath(self.qcbefore, rs['leftpretrim'])
                    rs['fastqcbeforeright'] = PipelineUtils.fastqcPath(self.qcbefore, rs['rightpretrim'])
                    if not self.dry:
                        ACT.submit("fastqc.qsub {} {}".format(rs['leftpretrim'], self.qcbefore), done="fqc.@.done")
                        ACT.submit("fastqc.qsub {} {}".format(rs['rightpretrim'], self.qcbefore), done="fqc.@.done")
                        nqc += 2

                if self.trim == 'trimmomatic':
                    LOG.log("Calling trimmomatic on `{}' and `{}'", rs['leftpretrim'], rs['rightpretrim'])
                    LOG.log("Trimmomatic main output to `{}' and `{}'", rs['left'], rs['right'])
                    if not self.dry:
                        tr = ACT.submit("trimmomatic.qsub {} {} {} {} {} {} {} {}".format(
                                rs["leftpretrim"], rs["rightpretrim"],
                                rs["left"], rs["trunpair1"], 
                                rs["right"], rs["trunpair2"],
                                adapt, minlen), done="trim.@.done")
                        ntrim += 1
                elif props['prog'] == 'sickle':
                    print "*** not implemented yet"

                if self.qcafter:
                    LOG.log("Performing QC on `{}' and `{}'", rs['left'], rs['right'])
                    rs['fastqcafterleft'] = PipelineUtils.fastqcPath(self.qcafter, rs['left'])
                    rs['fastqcafterright'] = PipelineUtils.fastqcPath(self.qcafter, rs['right'])
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
                    rs['fastqcbeforeleft'] = PipelineUtils.fastqcPath(self.qcbefore, rs['leftpretrim'])
                    if not self.dry:
                        ACT.submit("fastqc.qsub {} {}".format(rs['leftpretrim'], self.qcbefore), done="fqc.@.done")
                        nqc += 1

                if self.trim == 'trimmomatic':
                    LOG.log("Calling trimmomatic on `{}'", rs['leftpretrim'])
                    LOG.log("Trimmomatic main output to `{}'", rs['left'])
                    if not self.dry:
                        tr = ACT.submit("trimmomatic-se.qsub {} {} {} {}".format(rs["leftpretrim"], rs["left"], adapt, minlen),
                                        done="trim.@.done")
                        ntrim += 1
                elif props['prog'] == 'sickle':
                    print "*** not implemented yet"

                if self.qcafter:
                    LOG.log("Performing QC on `{}'", rs['left'])
                    rs['fastqcafterleft'] = PipelineUtils.fastqcPath(self.qcafter, rs['left'])
                    if not self.dry:
                        ACT.submit("fastqc.qsub {} {}".format(rs["left"], self.qcafter), after=tr, done="fqc.@.done")
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
        tbl2 = Table.ScrollingTable(id='tbl2', align="LRLRL", caption='Number of reads in input files and links to QC reports.', filename="table1.csv")
        tbl2.startHead()
        tbl2.addHeaderRow(["Sample", "Readset", "Reads before trim", "QC before trim", "Reads after trim", "QC after trim"])
        tbl2.startBody()
        for smp in SC.samples:
            tbl2.addRowHeader(smp['name'], rowspan=len(smp['readsets']) + 1)
            for rs in smp['readsets']:
                if rs['paired']:
                    tbl2.addRow([ Table.C(rs['name']), 
                                  Table.D(prereadcounts[rs['name']]),
                                  Table.C(rs['fastqcbeforeleft'] + "<br>" + rs['fastqcbeforeright']),
                                  Table.D(readcounts[rs['name']]),
                                  Table.C(rs['fastqcafterleft'] + "<br>" + rs['fastqcafterright'])
                                  ])
                else:
                    tbl2.addRow([ Table.C(rs['name']), 
                                  Table.D(prereadcounts[rs['name']]),
                                  Table.C(rs['fastqcbeforeleft']),
                                  Table.D(readcounts[rs['name']]),
                                  Table.C(rs['fastqcafterleft'])
                                  ])
            tbl2.addRow([ Table.C("<b>" + smp['name'] + "</b> (total)"),
                          Table.D(prereadcounts[smp['name']]),
                          Table.C("&nbsp;"),
                          Table.D(readcounts[smp['name']]), 
                          Table.P(1.0 * readcounts[smp['name']] / prereadcounts[smp['name']])
                          ])
        tbl2.toHTML(ACT.out)
        return True

### Bowtie2

class Bowtie2(Line):
    name = "Bowtie2"
    version = ""
    extra = ""                  # Extra bowtie2 options

    def Verify(self):
        #if self.dry:
        #    return True         # Don't stop the pipeline if we're not doing this
        ACT = self.actor
        verstring = ACT.shell("module load bowtie2; bowtie2 --version|head -1")
        p = verstring.find("version")
        if p > 0:
            self.version = verstring[p+8:]
            return True
        else:
            self.status = "Bowtie2 not found."
            return False

    def Setup(self):
        ACT = self.actor
        SC = ACT.sc

        self.extra = Utils.dget('extra', self.properties, self.extra)

        good = True
        ## Get name of Bowtie2 index for each sample
        for s in SC.samples:
            s['bt2idx'] = ACT.getConf("bt2idx", s['name']) or ACT.bt2idx
            if s['bt2idx'] == None:
                good = False

        if good or self.dry:
            return True
        else:
            self.status = "Bowtie2 index not specified. Please add a `bt2idx' entry in the configuration file."
            return False

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nbowtie = 0
        for smp in SC.samples:
            paired = False
            smp["bam"] = smp["name"] + ".bam"
            smp["bt2stats"] = smp["name"] + ".bt2stats.txt"
            lfq = []
            rfq = []
            for rs in smp['readsets']:
                lfq.append(rs['left'])
                if rs['paired']:
                    paired = True
                    rfq.append(rs['right'])
            lfq = ",".join(lfq)
            rfq = ",".join(rfq)
            LOG.log("Starting bowtie2 for sample {}, output to {}.", smp["name"], smp["bam"])
            if not self.dry:
                if paired:
                    cmdline = "bowtie2-pe.qsub sm={} stats={} {} {} {} {} {}".format(smp["name"], smp["bt2stats"], self.extra, smp["bt2idx"], lfq, rfq, smp["bam"])
                else:
                    cmdline = "bowtie2.qsub sm={} stats={} {} {} {} {}".format(smp["name"], smp["bt2stats"], self.extra, smp["bt2idx"], lfq, smp["bam"])
                ACT.submit(cmdline, done="bowtie2.@.done")
                nbowtie += 1
        return ACT.wait(("bowtie2.@.done", nbowtie))

    def Report(self):
        #print "printing bowtie report to {}".format(ACT.out)
        #raw_input()
        ACT = self.actor
        SC = ACT.sc

        readcounts = Utils.fileToDict(ACT.fastqCounts)

        ACT.scene("Alignment to genome")
        ACT.reportf("""The input sequences were aligned to the genome using Bowtie <b>{}</b>. The following table reports the number of aligned reads for each 
sample. The WIG files can be uploaded to the <A href='http://genome.ucsc.edu/'>UCSC
Genome Browser</A> as <A href='http://genome.ucsc.edu/cgi-bin/hgCustom'>custom tracks</A>.""", self.version)
        #geAlignments = Utils.fileToDict(ACT.afterCounts)
        tbl3 = Table.ScrollingTable(id='tbl3', align="RRC", caption='Number of alignments to genome.')
        tbl3.startHead()
        tbl3.addHeaderRow(["Sample", "Total reads", "Aligned reads", "Concordant alignment rate", "Bowtie2 report"])
        tbl3.startBody()
        for smp in SC.samples:
            sn = smp['name']
            fm = PipelineUtils.FileMatcherBT2()
            bt2 = fm.matchFile(smp["bt2stats"])
            tbl3.addRowHeader(sn)
            tbl3.addRow([ Table.D(bt2[1].groups()[0]),
                          Table.D(bt2[2].groups()[0]),
                          #Table.F(1.0 * Utils.dget(sn, geAlignments, default=0) / readcounts[sn]),
                          Table.P(bt2[2].groups()[1]),
                          Table.A(PipelineUtils.htmlize(smp['bt2stats'])) ])
        tbl3.toHTML(ACT.out)
        return True

### TopHat

class TopHat(Line):
    name = "TopHat"
    version = ""

    def Setup(self):
        ACT = self.actor
        SC = ACT.sc

        ## Get name of Bowtie2 index for each sample
        for s in SC.samples:
            s['bt2idx'] = ACT.getConf("bt2idx", s['name']) or ACT.bt2idx
        return True

    def Verify(self):
        ACT = self.actor
        version = ACT.shell("module load tophat; tophat -v")
        if version.startswith("TopHat"):
            self.version = version
            return True
        else:
            return False

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        ntophat = 0
        for smp in SC.samples:
            # Create an output directory for each sample
            outdir = smp['name'] + ".tophat"
            smp["bam"] = outdir + "/accepted_hits.bam" # Name of BAM file produced by TopHat
            ACT.mkdir(outdir)
            fqlistfile = self.tempfile(outdir + "/" + Utils.id_generator(10, prefix='tmp-'))
            lfq = []
            rfq = []
            for rs in smp['readsets']:
                lfq.append(rs['left'])
                if rs['paired']:
                    rfq.append(rs['right'])
            lfq = ",".join(lfq)
            rfq = ",".join(rfq)
            with open(fqlistfile, "w") as out:
                out.write("export LEFT={}\nexport RIGHT={}\n".format(lfq, rfq))
            if not self.dry:
                ACT.submit("tophat.qsub {} {} {}".format(smp["bt2idx"], fqlistfile, outdir), done="tophat.@.done")
                ntophat += 1
        return ACT.wait(("tophat.@.done", ntophat))

    def PostExecute(self):
        self.cleanTempfiles()

class Cufflinks(Line):
    """This object calls cufflinks."""
    name = "Cufflinks"
    version = ""

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        ncuff = 0
        with open("GTFS", "w") as out:
            for smp in SC.samples:
                # Reuse TopHat's output directory?
                outdir = smp['name'] + ".tophat"
                out.write(outdir + "/transcripts.gtf\n")
                if not self.dry:
                    ACT.submit("cuff.qsub cufflinks {} {} {}".format(outdir, smp['bam'], ACT.gtf), done="cuff.@.done")
                    ncuff += 1
        ACT.wait(("cuff.@.done", ncuff))
        ACT.submit("cuff.qsub cuffmerge merged_asm/ GTFS gtf={} seq={}".format(ACT.gtf, ACT.reference), done="cuffmerge.done")
        return ACT.wait("cuffmerge.done")

class Cuffdiff(Line):
    """This object calls cuffdiff."""
    name = "Cuffdiff"
    version = ""

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        ncdiff = 0
        for contrast in SC.contrasts:
            ctrlcond = contrast['control']
            ctrlSamples = SC.conditionSamples(ctrlcond)
            testcond = contrast['test']
            testSamples = SC.conditionSamples(testcond)
            outdir = contrast['cuffdiffdir'] = "{}.vs.{}.cuffdiff".format(testcond, ctrlcond)
            ACT.mkdir(outdir)
            if not self.dry:
                cmdline = "cuff.qsub cuffdiff {} {} {},{} {} {} fdr={}".format(outdir, "merged_asm/merged.gtf", ctrlcond, testcond,
                                                                               ",".join([ s['bam'] for s in ctrlSamples ]),
                                                                               ",".join([ s['bam'] for s in testSamples ]),
                                                                               ACT.fdr)
                ACT.submit(cmdline, done="cdiff.@.done")
                ncdiff += 1
        return ACT.wait(("cdiff.@.done", ncdiff))

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

            outdir = smp["name"]
            ACT.mkdir(outdir)
            smp["bam"] = smp["name"] + "/Aligned.sortedByCoord.out.bam" # Name of genome BAM file produced by STAR
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
                ACT.submit("star.qsub {} {} {} {}".format(smp["staridx"], lfq, rfq, outdir), done="star.@.done")
                nstar += 1
         
            # for rs in smp['readsets']:
            #     outdir = rs["name"]
            #     ACT.mkdir(outdir)
            #     rs["bam"] = rs["name"] + "/Aligned.sortedByCoord.out.bam" # Name of BAM file produced by STAR                                                               
            #     if not self.dry:
            #         if rs['paired']:
            #             LOG.log("Running STAR aligner on `{}' and `{}'", rs["left"], rs["right"])
            #             ACT.submit("star.qsub {} {} {} {}".format(smp["staridx"], rs["left"], rs["right"], outdir), done="star.@.done")
            #         else:
            #             LOG.log("Running STAR aligner on `{}'", rs["left"])
            #             ACT.submit("star.qsub {} {} {}".format(smp["staridx"], rs["left"], outdir), done="star.@.done")
            #         nstar += 1
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
        tbl3.addHeaderRow(["Sample", "Genome alignments", "Genome alignment rate", "Transcriptome alignments", "Transcriptome alignment rate", "Alignment report"])
        tbl3.startBody()
        for smp in SC.samples:
            logfile = smp['name'] + "/Log.final.out"
            loghtml = smp['name'] + "/Log.final.html"
            PipelineUtils.htmlize(logfile, loghtml)
            sn = smp['name']
            tbl3.addRowHeader(sn)
            tbl3.addRow([ Table.D(Utils.dget(sn, geAlignments, default=0)), 
                          Table.F(1.0 * Utils.dget(sn, geAlignments, default=0) / readcounts[sn]),
                          Table.D(txAlignments[sn]),
                          Table.F(1.0 * txAlignments[sn] / readcounts[sn]),
                          Table.A(loghtml, text=logfile)
                          ])
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

    def combineCounts(self, SC, log=None):
        idx = 1
        good = True
        out = open(self.outfile, "w") if self.outfile else None
        # self.actor.shell("ls -l {}*".format(prefix))
        try:
            for smp in SC.samples:
                smptotal = 0
                print "{} => {} readsets".format(smp['name'], len(smp['readsets']))
                for rs in smp['readsets']:
                    cntfile = "{}-{}.out".format(self.propname, idx)
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
                            if log:
                                pass # log.log("{}: {} reads", rs['name'], c)  ## Too verbose...
                    else:
                        self.status = "Warning: empty FASTQ file {}".format(cntfile)
                        good = False
                    idx += 1
                if self.key:
                    smp[self.propname] = smptotal
                if out:
                    out.write("{}\t{}\n".format(smp['name'], smptotal))
                if log:
                    log.log("{}: {} total reads", smp['name'], smptotal)
        finally:
            if out:
                out.close()
        return good

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        id = Utils.id_generator(10, prefix='tmp-')
        self.prefix = id
        fastqnames = self.tempfile(id + "-names.txt")
        qsub = self.tempfile(id + ".qsub")
        self.done = id + ".@.done"

        # Store the names of all FASTQ files in a file
        with open(fastqnames, "w") as out:
            for smp in SC.samples:
                for rs in smp['readsets']:
                    out.write("{}\n".format(rs['left']))
                    self.nfq += 1
        
        LOG.log("Submitting array job to compute number of reads in {} FASTQ files.", self.nfq)
        # Write the qsub script
        self.writeQsub(qsub, fastqnames, self.propname)
        if not self.dry:
            for i in range(1, self.nfq + 1):
                ACT.submit("{} {}".format(qsub, i), done=self.done)
        if not self.delay:
            ACT.wait((self.done, self.nfq))
            if self.outfile:
                LOG.log("Collecting FASTQ counts into {}.", self.outfile)
            if self.propname:
                LOG.log("Storing FASTQ counts in key {}.", self.propname)
            return self.combineCounts(SC, LOG) # Signal error if one or more FASTQ files are empty
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
        good = self.combineCounts(SC, LOG) # Signal error if one or more FASTQ files are empty

        #LOG.log("Deleting temporary count files.")
        #self.cleanTempfiles()
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
    propname = "bamcount"
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

    def combineCounts(self, SC, log=None):
        idx = 1
        good = True
        out = open(self.outfile, "w") if self.outfile else None
        # self.actor.shell("ls -l {}*".format(prefix))
        try:
            for smp in SC.samples:
                cntfile = "{}-{}.out".format(self.propname, idx)
                c = Utils.safeReadIntFromFile(cntfile, default=0, delay=2, maxtries=30)
                if c == 0:
                    self.status = "Warning: empty BAM file count in {}".format(cntfile)
                    good = False
                else:
                    if self.key:
                        smp[self.key] = c
                    if out:
                        out.write("{}\t{}\n".format(smp['name'], c))
                    if log:
                        log.log("{}: {} reads.", smp['name'], c)
                    idx += 1
        finally:
            if out:
                out.close()
        return good

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        id = Utils.id_generator(10, prefix='tmp-')
        self.prefix = id
        bamnames = self.tempfile(id + "-names.txt")
        qsub = id + ".qsub"
        self.done = id + ".@.done"

        # Store the names of all BAM files in a file
        with open(bamnames, "w") as out:
            for smp in SC.samples:
                out.write("{}\n".format(smp[self.source]))
                self.nbams += 1
        
        LOG.log("Submitting array job to compute number of reads in {} BAM files.", self.nbams)
        # Write the qsub array script
        self.writeQsub(qsub, bamnames, self.propname)
        if not self.dry:
            for i in range(1, self.nbams+1):
                ACT.submit("{} {}".format(qsub, i), done=self.done)
        if self.delay:
            return True         # We'll collect results in PostExecute
        else:
            ACT.wait((self.done, self.nbams))
            if self.outfile:
                LOG.log("Collecting BAM counts into {}.", self.outfile)
            if self.propname:
                LOG.log("Storing BAM counts in key {}.", self.propname)
            return self.combineCounts(SC, LOG) # Signal error if one or more BAM files are empty

    def PostExecute(self):
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log

        good = True
        if self.delay:
            if not self.dry:
                ACT.wait((self.done, self.nbams))
            if self.outfile:
                LOG.log("Collecting BAM counts into {}.", self.outfile)
            if self.propname:
                LOG.log("Storing BAM counts in key {}.", self.propname)
            good = self.combineCounts(SC, LOG) # Signal error if one or more BAM files are empty
        
        LOG.log("Deleting temporary count files.")
        self.cleanTempfiles()
        return good

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
        if not self.dry:
            for smp in SC.samples:
                # Utils.safeRemoveFile(smp['origbam'])
                Utils.safeRemoveFile(smp['markdup'])
        return True

### BAM merger 

class BAMmerger(Line):
    """Merge sets of BAM files with bamtools."""

    name = "BAM merger"
    byCondition = False
    indexBAM = False
    removeOriginal = False
    bamkey = 'bam'

    def Setup(self):
        self.byCondition = Utils.dget('byCondition', self.properties, self.byCondition)
        self.indexBAM = Utils.dget('indexBAM', self.properties, self.indexBAM)
        self.removeOriginal = Utils.dget('remove', self.properties, self.removeOriginal)
        self.bamkey = Utils.dget('bamkey', self.properties, self.bamkey)
        return True

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nmerge = 0
        nindex = 0
        for smp in SC.samples:
            newbamname = smp['name'] + ".bam"
            if self.bamkey in smp:    # We already have a per-sample BAM file
                if smp[self.bamkey] != newbamname: # but it doesn't have the name we want...
                    LOG.log("Copying {} to {}", smp[self.bamkey], newbamname)
                    if ACT.missingOrStale(newbamname, other=smp[self.bamkey]):
                        ACT.copy(smp[self.bamkey], newbamname)
                    smp[self.bamkey] = newbamname
                    if self.indexBAM and not self.dry and ACT.missingOrStale(smp[self.bamkey] + ".bai", other=smp[self.bamkey]):
                        ACT.submit("bam-index.qsub {}".format(smp[self.bamkey]), done='index.@.done')
                        nindex += 1
            else:
                LOG.log("Merging BAM files for sample {} into {}", smp['name'], smp[self.bamkey])
                cmdline = "generic.qsub bamtools merge -out " + smp[self.bamkey]
                needed = False
                for rs in smp['readsets']:
                    if 'bam' in rs:
                        needed = True
                        cmdline += " -in " + rs[self.bamkey]
                if needed:
                    cmdline += " module:bamtools"
                    if not self.dry:
                        j = ACT.submit(cmdline, done='merge.@.done')
                        nmerge += 1
                        if self.indexBAM and ACT.missingOrStale(smp[self.bamkey] + ".bai", other=smp[self.bamkey]):
                            ACT.submit("bam-index.qsub {}".format(smp[self.bamkey]), after=j, done='index.@.done')
                            nindex += 1
                ACT.wait(('merge.@.done', nmerge))

        if self.byCondition:
            nmerge = 0
            nindex = 0
            for cond in SC.conditions:
                cond['bam'] = cond['name'] + ".bam"
                LOG.log("Merging BAM files for condition {} into {}", cond['name'], cond['bam'])
                j = False
                condbams = SC.conditionBAMs(cond['name'], key=self.bamkey)
                if ACT.missingOrStale(cond['bam'], condbams):
                    cmdline = "generic.qsub bamtools merge -out " + cond['bam']
                    for sb in condbams:
                        cmdline += " -in " + sb
                    cmdline += " module:bamtools"
                    if not self.dry:
                        ACT.delete(cond['bam'] + ".bai") # If we're merging, we need to re-index for sure
                        j = ACT.submit(cmdline, done='cmerge.@.done')
                        nmerge += 1
                if not self.dry and self.indexBAM:
                    if j:
                        ACT.submit("bam-index.qsub {}".format(cond['bam']), after=j, done='index.@.done')
                        nindex += 1
                    elif ACT.missingOrStale(cond['bam'] + ".bai", other=cond['bam']):
                        ACT.submit("bam-index.qsub {}".format(cond['bam']), done='index.@.done')
                        nindex += 1

            return ACT.wait(('cmerge.@.done', nmerge), ('index.@.done', nindex))
        else:
            return True

### BAM concatenator

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

### BAM coverage

class BAMcoverage(Line):
    """Compute BAM coverage statistics."""
    name = "BAM coverage"
    ncovs = 0                   # Number of bamcoverage processes running

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        self.ncovs = 0
        for smp in SC.samples:
            if 'genomebam' in smp:
                key = 'genomebam'
            elif 'bam' in smp:
                key = 'bam'
            else:
                continue
            cov = smp['bamcov'] = smp['name'] + ".cov.csv"
            if ACT.missingOrStale(cov, other=smp[key]) and not self.dry:
                ACT.submit("bamcoverage.qsub {} {} Y {}".format(smp[key], cov, ACT.bamcovd), done="bamcov.@.done")
                self.ncovs += 1
        for cond in SC.conditions:
            if 'bam' in cond:
                cov = cond['bamcov'] = cond['name'] + ".cov.csv"
                if ACT.missingOrStale(cov, other=cond['bam']) and not self.dry:
                    ACT.submit("bamcoverage.qsub {} {} Y {}".format(cond['bam'], cov, ACT.bamcovd), done="bamcov.@.done")
                    self.ncovs += 1
        return True

    def PostExecute(self):
        ACT = self.actor
        SC = ACT.sc

        ACT.wait(("bamcov.@.done", self.ncovs))
        ACT.samplecov = "{}.sample.cov.xlsx".format(ACT.Name)
        cmdline = "module load dibig_tools; csvtoxls.py " + ACT.samplecov
        for smp in SC.samples:
            cmdline += " {} -name {} -firstrowhdr".format(smp['bamcov'], smp['name'])
        ACT.shell(cmdline)

        ACT.condcov = "{}.cond.cov.xlsx".format(ACT.Name)
        cmdline = "module load dibig_tools; csvtoxls.py " + ACT.condcov
        for cond in SC.conditions:
            if 'bamcov' in cond:
                cmdline += " {} -name {} -firstrowhdr".format(cond['bamcov'], cond['name'])
        ACT.shell(cmdline)
        return True
    
    def Report(self):
        ACT = self.actor
        SC = ACT.sc
        
        ACT.scene("Genome coverage")
        ACT.reportf("""The following table reports the overall and effective genome coverage in each sample. The <i>Total nt</i> column reports
the total number of nucleotides sequenced, i.e. the number of aligned reads times the length of each read. <i>Coverage</i> is this number 
divided by the size of the genome. <i>Effective bp</i> reports the number of bases in the genome having coverage greater than <b>{}</b>, and the
<i>Effective Perc</i> column shows what percentage this is of the genome size. Note that, especially in the case of RNA-seq, the effective
genome size may be much smaller than the full size. <i>Eff Coverage</i> is the average coverage over the effectively covered fraction of
the genome. """.format(ACT.bamcovd))
        tbl = Table.ScrollingTable(id='tblscov', caption="Genome coverage by sample.")
        tbl.startHead()
        tbl.addHeaderRow(['Name', 'Total nt', 'Coverage', 'Effective bp', 'Effective Perc', 'Eff Coverage'])
        tbl.startBody()
        for smp in SC.samples:
            if 'bamcov' in smp:
                data = Utils.fileToList(smp['bamcov'])[0]
                row = data[-1]  #  we want the last entry (total)
                tbl.addRow([ Table.C(smp['name']), 
                             Table.D(int(row[1])), 
                             Table.F(float(row[3])), 
                             Table.D(int(row[4])), 
                             Table.P(float(row[5]) / 100.0), 
                             Table.F(float(row[6])) ])
        tbl.toHTML(ACT.out)

        ACT.reportf("""The following table reports the overall and effective genome coverage in each condition.""")
        tbl = Table.ScrollingTable(id='tblccov', caption="Genome coverage by condition")
        tbl.startHead()
        tbl.addHeaderRow(['Name', 'Total nt', 'Coverage', 'Effective bp', 'Effective Perc', 'Eff Coverage'])
        tbl.startBody()

        for cond in SC.conditions:
            if 'bamcov' in cond:
                data = Utils.fileToList(cond['bamcov'])[0]
                row = data[-1]  #  we want the last entry (total)
                tbl.addRow([ Table.C(cond['name']),
                             Table.D(int(row[1])), 
                             Table.F(float(row[3])), 
                             Table.D(int(row[4])), 
                             Table.P(float(row[5]) / 100.0), 
                             Table.F(float(row[6])) ])
                
        tbl.toHTML(ACT.out)

        ACT.file(ACT.samplecov, description="Per-chromosome coverage data, by sample.")
        ACT.file(ACT.condcov, description="Per-chromosome coverage data, by condition.")
                 
        return True

### Kallisto

class Kallisto(Line):
    """Perform transcriptome quantification using Kallisto. Calls kallisto.qsub on the BAM files for each sample."""
    name = "Kallisto"
    version = ""
    index = False
    
    def Setup(self):
        ACT = self.actor
        self.index = ACT.getConf("kallistoidx") or self.index
        return self.index

    def Verify(self):
        ACT = self.actor
        version = ACT.shell("module load kallisto; kallisto version")
        if version.startswith("kallisto"):
            p = version.rfind(" ")
            self.version = version[p+1:]
            return True
        else:
            return False

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nkall = 0
        for smp in SC.samples:
            smp['kallistodir'] = outdir = smp['name'] + ".kallisto.d"
            ACT.mkdir(outdir)
            cmdline = "kallisto.qsub {} {}".format(self.index, outdir)
            for rs in smp['readsets']:
                cmdline += " " + rs['left'] + " " + rs['right']
            if not self.dry:
                ACT.submit(cmdline, done="kall.@.done")
                nkall += 1
        ACT.wait(("kall.@.done", nkall))
        nmat = 0
        if not self.dry:
            ACT.submit(writeMatrixScript("allgenematrix.qsub", SC.samples, "genes.rawmatrix.in.csv", mode='k'), done="raw.@.done")
            nmat += 1
        ACT.wait(("raw.@.done", nmat))
        if not self.dry:
            PipelineUtils.annotateFile("genes.rawmatrix.in.csv", "genes.rawmatrix.csv", ACT.txtable, annot=['gene_name', 'gene_biotype'], annotNames=['Gene', 'Biotype'])
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
        ACT.wait(("rsem.@.done", nrsem))

        nmat = 0
        if  not self.dry:
            ACT.submit(writeMatrixScript("allgenematrix.qsub", SC.samples, "genes.rawmatrix.in.csv", mode='g', column='FPKM'), done="raw.@.done")
            ACT.submit(writeMatrixScript("allisomatrix.qsub", SC.samples, "transcripts.rawmatrix.in.csv", mode='i', column='FPKM'), done="raw.@.done")
            nmat = 2
        ACT.wait(("raw.@.done", nmat))
        if not self.dry:
            PipelineUtils.annotateFile("genes.rawmatrix.in.csv", "genes.rawmatrix.csv", 
                                       ACT.gtable, annot=['gene_name', 'gene_biotype'], annotNames=['Gene', 'Biotype'])
            PipelineUtils.annotateFile("transcripts.rawmatrix.in.csv", "transcripts.rawmatrix.csv",
                                       ACT.txtable, annot=['transcript_name', 'gene_name', 'gene_biotype'], annotNames=['Transcript', 'Gene', 'Biotype'])
        return True

    def exprScatterplots(self, ACT):
        SC = ACT.sc

        # *** Hack: get sample names from file header
        with open("allExpressions.csv", "r") as f:
            hdr = f.readline().rstrip("\r\n").split("\t")
        for i in range(len(hdr)):
            if hdr[i].endswith(".genes.results"):
                hdr[i] = hdr[i][:-14]

        for cond in SC.conditions:
            samples = SC.conditionSamples(cond)
            # print samples
            s0 = samples[0]
            s1 = samples[1]
            idx0 = hdr.index(s0['name'])
            idx1 = hdr.index(s1['name'])
            # print "{}=>{}, {}=>{}".format(s0['name'], idx0, s1['name'], idx1)
            cond['scatterplot'] = cond['name'] + ".scatt.png"
            p = exprScatterplot(cond['scatterplot'], "allExpressions.csv", [idx0], [idx1], label1=s0['name'], label2=s1['name'], fc=None)
            cond['correlation'] = p.correlation[0][1]

    def showScatterplots(self, ACT):
        SC = ACT.sc
        out = ACT.out
        out.write("""<CENTER><TABLE class='figure'><TR><TD class='figure'>\n""")
        out.write("""  <TABLE width='100%'><TR><TD align='left' valign='top'>\n""")
        out.write("""    <TABLE class='legend'><TH class='legend'>Condition</TH><TH class='legend'>Correlation</TH><TH class='legend'>Corr Score</TH></TR>\n""")
        img1 = SC.conditions[0]['scatterplot']
        for cond in SC.conditions:
            score = -math.log10(1-cond['correlation'])
            out.write("""    <TR onclick="document.getElementById('condscatter').src='{}';"><TD class='legend'>{}</TD><TD align='right' class='legend'>{}</TD><TD align='right' class='legend'>{}</TD></TR>\n""".format(cond['scatterplot'],
                    cond['name'], cond['correlation'], Utils.f3dd(score)))
        out.write("""    </TABLE>
  </TD><TD align='right' valign='top'>
    <A href='{}'><IMG class='bordered' id='condscatter' src='{}' width='500'/></A>
  </TD></TR></TABLE>
</TD></TR>
</TABLE></CENTER>
""".format(img1, img1))

    def Report(self):
        ACT = self.actor
        SC = ACT.sc

        self.exprScatterplots(ACT)

        ACT.scene("Expression analysis - quantification")
        ACT.reportf("""Gene and transcript expression values were quantified using <b>{}</b>. The following files contain the raw FPKM values for all
genes/transcripts in all samples. <b>NOTE:</b> these values are not normalized yet, please apply the appropriate normalization before using them in
analysis.""".format(self.version))
        ACT.file("genes.rawmatrix.csv", description="Matrix of FPKM values for all genes in all samples.")
        ACT.file("transcripts.rawmatrix.csv", description="Matrix of FPKM values for all transcripts in all samples.")

        self.showScatterplots(ACT)
        return True

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
        contrast['grsemqsub'] = gfilename = "{}.vs.{}.g.qsub".format(testcond, ctrlcond)
        contrast['gmatrix'] = "{}.vs.{}.gmatrix.csv".format(testcond, ctrlcond)
        contrast['gdiff'] = "{}.vs.{}.gdiff.csv".format(testcond, ctrlcond)
        contrast['gfdr'] = "{}.vs.{}.gfdr.csv".format(testcond, ctrlcond)

        writeContrastMatrix(contrast, testSamples + ctrlSamples, ntest, nctrl, ACT.fdr, erccdb=ACT.erccdb, mode="g")

        contrast['irsemqsub'] = ifilename = "{}.vs.{}.i.qsub".format(testcond, ctrlcond)
        contrast['imatrix'] = "{}.vs.{}.imatrix.csv".format(testcond, ctrlcond)
        contrast['idiff'] = "{}.vs.{}.idiff.csv".format(testcond, ctrlcond)
        contrast['ifdr'] = "{}.vs.{}.ifdr.csv".format(testcond, ctrlcond)

        writeContrastMatrix(contrast, testSamples + ctrlSamples, ntest, nctrl, ACT.fdr, erccdb=ACT.erccdb, mode="i")

        return (gfilename, ifilename)

    def doExprScatterplots(self, ACT, contrast):
        SC = ACT.sc
        ctrlcond = contrast['control']
        ctrlSamples = SC.conditionSamples(ctrlcond)
        testcond = contrast['test']
        testSamples = SC.conditionSamples(testcond)
        nctrl = len(ctrlSamples)
        ntest = len(testSamples)
        ctrlidx = [x+1       for x in range(nctrl)]
        testidx = [x+nctrl+1 for x in range(ntest)]
        contrast['scatterplot'] = "{}.vs.{}.png".format(testcond, ctrlcond)
        # print (contrast['scatterplot'], contrast['gmatrix'], ctrlidx, testidx, ctrlcond, testcond)
        p = exprScatterplot(contrast['scatterplot'], contrast['gmatrix'], testidx, ctrlidx, label1=testcond, label2=ctrlcond, fc=ACT.fc)
        contrast['correlation'] = p.correlation[0][1]

    def showScatterplots(self, ACT):
        SC = ACT.sc
        out = ACT.out
        out.write("""<CENTER><TABLE class='figure'><TR><TD class='figure'>\n""")
        out.write("""  <TABLE width='100%'><TR><TD align='left' valign='top'>\n""")
        out.write("""    <TABLE class='legend'><TH class='legend'>Test</TH><TH class='legend'>Ctrl</TH><TH class='legend'>Correlation</TH><TH class='legend'>Corr Score</TH></TR>\n""")
        img1 = SC.contrasts[0]['scatterplot']
        for contrast in SC.contrasts:
            score = -math.log10(1-contrast['correlation'])
            out.write("""    <TR onclick="document.getElementById('scatter').src='{}';"><TD class='legend'>{}</TD><TD class='legend'>{}</TD><TD align='right' class='legend'>{}</TD><TD align='right' class='legend'>{}</TD></TR>\n""".format(contrast['scatterplot'],
                    contrast['test'], contrast['control'], contrast['correlation'], Utils.f3dd(score)))
        out.write("""    </TABLE>
  </TD><TD align='right' valign='top'>
    <A href='{}'><IMG class='bordered' id='scatter' src='{}' width='500'/></A>
  </TD></TR></TABLE>
</TD></TR>
</TABLE></CENTER>
""".format(img1, img1))

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

    def PostExecute(self):
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log

        allfinalg = []
        allfinalc = []
        allfinali = []
        all2finalg = []
        all2finalc = []
        all2finali = []
        mergecmdline = "module load dibig_tools; rnaseqtools.py merge "
        allexpcmdline = "module load dibig_tools; rnaseqtools.py allexp -o allExpressions.in.csv "
        LOG.log("Generating final differential expression tables.")
        for contr in SC.contrasts:
            self.doExprScatterplots(ACT, contr)
            allexpcmdline += " " + contr['gmatrix']
            ctrlcond = contr['control']
            testcond = contr['test']
            vs = "{}.vs.{}".format(testcond, ctrlcond)
            contr['genefinal']      = vs + ".geneDiff.csv"
            contr['codingfinal']    = vs + ".codinggeneDiff.csv"
            contr['isofinal']       = vs + ".isoDiff.csv"
            contr['allgenefinal']   = vs + ".allGeneDiff.csv"
            contr['allcodingfinal'] = vs + ".allCodingDiff.csv"
            contr['allisofinal']    = vs + ".allIsoDiff.csv"
            mergecmdline += vs + " "
            allfinalg.append(contr['genefinal'])
            allfinalc.append(contr['codingfinal'])
            allfinali.append(contr['isofinal'])
            all2finalg.append(contr['allgenefinal'])
            all2finalc.append(contr['allcodingfinal'])
            all2finali.append(contr['allisofinal'])
            contr['ngdiff'] = filterDiff(contr['gfdr'], contr['genefinal'], ACT.fc, translation=ACT.gtable,
                                         wanted=['gene_id', 'gene_biotype', 'gene_name'],
                                         wantedNames=["ENSG", "Biotype", "Gene"])
            contr['ncdiff'] = filterDiff(contr['gfdr'], contr['codingfinal'], ACT.fc, translation=ACT.gtable,
                                         wanted=['gene_id', 'gene_biotype', 'gene_name'],
                                         wantedNames=["ENSG", "Biotype", "Gene"],
                                         biotype='protein_coding')
            contr['nidiff'] = filterDiff(contr['ifdr'], contr['isofinal'], ACT.fc, translation=ACT.txtable,
                                         wanted=['transcript_id', 'gene_biotype', 'transcript_name', 'gene_name'],
                                         wantedNames=["ENST", "Biotype", "Transcript", "Gene"])
            contr['agdiff'] = filterDiff(contr['gdiff'], contr['allgenefinal'], 0, translation=ACT.gtable,
                                         wanted=['gene_id', 'gene_biotype', 'gene_name'],
                                         wantedNames=["ENSG", "Biotype", "Gene"])
            contr['acdiff'] = filterDiff(contr['gdiff'], contr['allcodingfinal'], 0, translation=ACT.gtable,
                                         wanted=['gene_id', 'gene_biotype', 'gene_name'],
                                         wantedNames=["ENSG", "Biotype", "Gene"],
                                         biotype='protein_coding')
            contr['aidiff'] = filterDiff(contr['gdiff'], contr['allisofinal'], 0, translation=ACT.gtable,
                                         wanted=['gene_id', 'gene_biotype', 'gene_name'],
                                         wantedNames=["ENSG", "Biotype", "Gene"])
        LOG.log("Executing: {}", mergecmdline)
        ACT.shell(mergecmdline)
        LOG.log("Executing: {}", allexpcmdline)
        ACT.shell(allexpcmdline)
        os.rename("merged.geneDiff.csv", "merged.geneDiff.in.csv")
        os.rename("merged.codinggeneDiff.csv", "merged.codinggeneDiff.in.csv")
        os.rename("merged.isoDiff.csv", "merged.isoDiff.in.csv")
        PipelineUtils.annotateFile("allExpressions.in.csv", "allExpressions.csv", ACT.gtable, annot=['gene_name', 'gene_biotype'], annotNames=['Gene', 'Biotype'])
        PipelineUtils.annotateFile("merged.geneDiff.in.csv", "merged.geneDiff.csv", ACT.gtable, annot=['gene_name', 'gene_biotype'], annotNames=['Gene', 'Biotype'])
        PipelineUtils.annotateFile("merged.codinggeneDiff.in.csv", "merged.codinggeneDiff.csv", ACT.gtable, annot=['gene_name', 'gene_biotype'], annotNames=['Gene', 'Biotype'])
        PipelineUtils.annotateFile("merged.isoDiff.in.csv", "merged.isoDiff.csv", ACT.txtable, annot=['transcript_name', 'gene_name', 'gene_biotype'], annotNames=['Transcript', 'Gene', 'Biotype'])

        ACT.ff = {}
        for v in ["allExpressions", "merged.allDiff", "genediff", "codingdiff", "isodiff", "allgenediff", "allcodingdiff", "allisodiff"]:
            ACT.ff[v] = ACT.Name + "-" + v + ".xlsx"

        with open("toxls.sh", "w") as out:
            out.write("#!/bin/bash\n\nmodule load dibig_tools\n\n")
            if ACT.missingOrStale(ACT.ff['allExpressions'], other="allExpressions.csv"):
                out.write("echo Creating {}...\n".format(ACT.ff['allExpressions']))
                out.write("csvtoxls.py {} -q allExpressions.csv\n".format(ACT.ff['allExpressions']))
            if ACT.missingOrStale(ACT.ff["merged.allDiff"], other=["merged.geneDiff.csv", "merged.codinggeneDiff.csv", "merged.isoDiff.csv"]):
                out.write("echo Creating {}...\n".format(ACT.ff['merged.allDiff']))
                out.write("csvtoxls.py {} -q merged.geneDiff.csv merged.codinggeneDiff.csv merged.isoDiff.csv\n".format(ACT.ff["merged.allDiff"]))
            if ACT.missingOrStale(ACT.ff["genediff"], other=allfinalg):
                out.write("echo Creating {}...\n".format(ACT.ff['genediff']))
                out.write("csvtoxls.py {} -q {}\n".format(ACT.ff["genediff"], " ".join(allfinalg)))
            if ACT.missingOrStale(ACT.ff["codingdiff"], other=allfinalc):
                out.write("echo Creating {}...\n".format(ACT.ff['codingdiff']))
                out.write("csvtoxls.py {} -q {}\n".format(ACT.ff["codingdiff"], " ".join(allfinalc)))
            if ACT.missingOrStale(ACT.ff["isodiff"], other=allfinali):
                out.write("echo Creating {}...\n".format(ACT.ff['isodiff']))
                out.write("csvtoxls.py {} -q {}\n".format(ACT.ff["isodiff"], " ".join(allfinali)))
            if ACT.missingOrStale(ACT.ff["allgenediff"], other=all2finalg):
                out.write("echo Creating {}...\n".format(ACT.ff['allgenediff']))
                out.write("csvtoxls.py {} -q {}\n".format(ACT.ff["allgenediff"], " ".join(all2finalg)))
            if ACT.missingOrStale(ACT.ff["allcodingdiff"], other=all2finalc):
                out.write("echo Creating {}...\n".format(ACT.ff['allcodingdiff']))
                out.write("csvtoxls.py {} -q {}\n".format(ACT.ff["allcodingdiff"], " ".join(all2finalc)))
            if ACT.missingOrStale(ACT.ff["allisodiff"], other=all2finali):
                out.write("echo Creating {}...\n".format(ACT.ff['allisodiff']))
                out.write("csvtoxls.py {} -q {}\n".format(ACT.ff["allisodiff"], " ".join(all2finali)))
        # ACT.message("*** skipping toxlsx.sh!")
        ACT.shell("source ./toxls.sh")

        # Compute intersections
        ACT.intersect = Intersector.Intersector(SC.contrasts, "codingfinal")
        ACT.intersect.splitUpDown()
        ACT.intersect.writeIntersectScript("int.sh")
        ACT.shell("source ./int.sh")
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
            tbl4.addRow([ Table.C(contr['test']), 
                          Table.C(contr['control']), 
                          Table.D(ng[1], cls='upreg'), 
                          Table.D(ng[2], cls='dnreg'),
                          Table.D(ng[0]),
                          Table.A(contr['codingfinal']),
                          Table.A(contr['gmatrix'])])
        tbl4.toHTML(ACT.out)
        ACT.file(ACT.ff["codingdiff"], description="Excel file containing differentially expressed genes for all contrasts (one sheet per contrast). Only includes coding genes.")
        ACT.file(ACT.ff["allcodingdiff"], description="Excel file containing differential expression values for all genes in all contrasts (one sheet per contrast). Only includes coding genes.")

        ACT.reportf("""The following table shows the number of common genes between each pair of contrasts.""")
        ACT.intersect.toHTML(ACT, 'tblint')

        self.showScatterplots(ACT)

        ACT.reportf("""The following table reports results from the same differential analysis as above, but includes all biotypes instead of coding genes only.""")

        tbl4b = Table.ScrollingTable(id='tbl4b', align="LLRRRCC", caption="Results of gene-level differential expression analysis (all biotypes).")
        tbl4b.startHead()
        tbl4b.addHeaderRow(["Test", "Control", "Overexpressed", "Underexpressed", "Total", "Table", "Expressions"])
        tbl4b.startBody()
        for contr in SC.contrasts:
            ng = contr['ngdiff']
            tbl4b.addRow([ Table.C(contr['test']), 
                           Table.C(contr['control']), 
                           Table.D(ng[1], cls='upreg'),
                           Table.D(ng[2], cls='dnreg'),
                           Table.D(ng[0]),
                           Table.A(contr['genefinal']),
                           Table.A(contr['gmatrix'])])
        tbl4b.toHTML(ACT.out)
        ACT.file(ACT.ff["genediff"], description="Excel file containing differentially expressed genes for all contrasts (one sheet per contrast). Includes all genes and pseudo-genes.")
        ACT.file(ACT.ff["allgenediff"], description="Excel file containing differential expression values for all genes in all contrasts (one sheet per contrast). Includes all genes and pseudo-genes.")
        ACT.file(ACT.ff["allExpressions"], description="Excel file containing normalized expression values for all genes in all conditions.")

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
            tbl5.addRow([ Table.C(contr['test']), 
                          Table.C(contr['control']), 
                          Table.D(ni[1], cls='upreg'),
                          Table.D(ni[2], cls='dnreg'),
                          Table.D(ni[0]),
                          Table.A(contr['isofinal']),
                          Table.A(contr['imatrix']) ])
        tbl5.toHTML(ACT.out)
        ACT.file(ACT.ff["isodiff"], description="Excel file containing differentially expressed isoforms for all contrasts (one sheet per contrast).")
        ACT.file(ACT.ff["allisodiff"], description="Excel file containing differential expression values for all isoforms in all contrasts (one sheet per contrast).")

        ACT.scene("Differential expression - combined files")
        ACT.reportf("""The following file contains <i>merged</i> differential expression data. The first sheet contains fold changes for all genes that were
found to be differentially expressed in at least one contrast. The second and third sheets contain the same information for coding genes only, and all
transcripts.""")
        ACT.file(ACT.ff["merged.allDiff"], description="Merged fold changes for all differentially expressed genes, coding genes, and transcripts respectively.")

        return True

### BAMtoWig

class BAMtoWig(Line):
    """Convert BAM files for samples to Wig files."""
    name = "BAM to Wig"
    window = 100
    scale = 1000000000
    countsfile = "afterCounts.csv"
    subtitle = "transcriptome"
    bamkey = "bam"

    def Setup(self):
        self.bamkey = Utils.dget('bamkey', self.properties, self.bamkey)
        self.scale = int(Utils.dget('scale', self.properties, self.scale))
        self.window = int(Utils.dget('window', self.properties, self.window))
        self.countsfile = Utils.dget('countsfile', self.properties, self.countsfile)
        self.subtitle = Utils.dget('subtitle', self.properties, self.subtitle)
        return True

    def Execute(self):
        #def PostExecute(self):      # done in postexecute so we'll have bam counts *** NO: replace with idxstats
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        #LOG.log("Loading BAM counts from {}", self.countsfile)
        #counts = Utils.fileToDict(self.countsfile)
        nwig = 0
        for smp in SC.samples:
            name = smp['name']
            smp['wigfile'] = name + ".wig"
            LOG.log("Writing WIG file for {} to {}", smp[self.bamkey], smp['wigfile'])
            if not self.dry:
                nreads = Utils.countReadsInBAM(smp[self.bamkey])
                LOG.log("Number of reads in {}: {}", smp[self.bamkey], nreads)
                ACT.submit("""generic.qsub bamToWig.py -o {} -n {} -w {} -s {} -t {} -d \"{}\" {} module:dibig_tools module:samtools""".format(
                        smp['wigfile'], nreads, self.window, self.scale, name, self.subtitle, smp[self.bamkey]), done="wig.@.done")
                nwig += 1

        # Do we also have per-condition BAM files?
        for cond in SC.conditions:
            if 'bam' in cond and os.path.isfile(cond['bam']):
                name = cond['name']
                cond['wigfile'] = name + ".wig"
                LOG.log("Writing WIG file for {} to {}", cond['bam'], cond['wigfile'])
                if not self.dry:
                    nreads = Utils.countReadsInBAM(cond['bam'])
                    LOG.log("Total number of reads for {}: {}", name, nreads)
                    ACT.submit("""generic.qsub bamToWig.py -o {} -n {} -w {} -s {} -t {} -d \"{}\" {} module:dibig_tools module:samtools""".format(
                            cond['wigfile'], nreads, self.window, self.scale, name, self.subtitle, cond['bam']), done="wig.@.done")
                    nwig += 1

        return ACT.wait(("wig.@.done", nwig))

### MACS

class MACScallpeak(Line):
    """Run MACS in callpeak mode."""
    name = "MACS (callpeak)"

    def Verify(self):
        ACT = self.actor
        verstring = ACT.shell("module load macs; macs2 --version 2>&1")
        if verstring[:5] == "macs2":
            self.version = verstring[6:]
            return True
        else:
            self.status = "MACS2 not found."
            return False

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
            contrast['macsdir'] = macsdir = "{}.vs.{}.macs".format(testcond, ctrlcond)
            cbams = ",".join(SC.conditionBAMs(ctrlcond))
            tbams = ",".join(SC.conditionBAMs(testcond))
            if ACT.fdr:
                pval = "pval={}".format(ACT.fdr)
            else:
                pval = ""
            if not self.dry:
                ACT.submit("macs-diff.qsub {} {} {} {} {}".format(contrast['name'], macsdir, tbams, cbams, pval), done="macs.@.done")
                nmacs += 1
        ACT.wait(("macs.@.done", nmacs))

        # Create a CSV version of the peaks file for frip calculation
        for cond in SC.conditions:
            # *** HACK: 25 is hardcoded, should be computed!
            ACT.shell("tail -n +25 {} > {}".format(cond['macsdir'] + "/" + cond['name'] + "_peaks.xls", cond['macsdir'] + "/" + cond['name'] + "_peaks.csv"))

        for contrast in SC.contrasts:
            ctrlcond = contrast['control']
            testcond = contrast['test']
            name = "{}.vs.{}".format(testcond, ctrlcond)
            pileup = contrast['macsdir'] + "/" + name + "_treat_pileup.bdg"
            pileup2 = contrast['macsdir'] + "/" + name + ".bedGraph"
            contrast['macsPileup'] = pileup2
            if os.path.isfile(pileup) and ACT.missingOrStale(pileup2, other=pileup):
                with open(pileup2, "w") as out:
                    out.write("track type=bedGraph\n")
                ACT.shell("psort -k1,1 -k2,2n {} >> {}".format(pileup, pileup2))
            summits = contrast['macsdir'] + "/" + name + "_summits.bed"
            summits2 = contrast['macsdir'] + "/" + name + ".summits.bedGraph"
            contrast['macsSummits'] = summits2
            if os.path.isfile(summits) and ACT.missingOrStale(summits2, other=summits):
                with open(summits2, "w") as out:
                    out.write("track type=bedGraph\n")
                ACT.shell("cut -f 1,2,3,5 {} | psort -k1,1 -k2,2n | mergeRegions.py >> {}".format(summits, summits2))
            narrow = contrast['macsdir'] + "/" + name + "_peaks.narrowPeak"
            narrow2 = contrast['macsdir'] + "/" + name + ".npeaks.bedGraph"
            contrast['macsNarrow'] = narrow2
            if os.path.isfile(narrow) and ACT.missingOrStale(narrow2, other=narrow):
                with open(narrow2, "w") as out:
                    out.write("track type=bedGraph\n")
                ACT.shell("cut -f 1,2,3,5 {} | psort -k1,1 -k2,2n | mergeRegions.py >> {}".format(narrow, narrow2))

        return True

    def Report(self):
        ACT = self.actor
        SC = ACT.sc

        ACT.scene("Peak detection")
        ACT.reportf("""Peak detection was performed using <B>MACS</B> version <b>{}</b>. The following table provides links to the <i>Pileup</i>, <i>narrowPeaks</i>, and <i>Summits</i> files for each contrast.
All files are in bedGraph format.""".format(self.version))
        tbl = Table.ScrollingTable(id='tblmacs', align="LLCCC", caption="Results of peak detection with MACS.")
        tbl.startHead()
        tbl.addHeaderRow(["Test", "Control", "Pileup", "Peaks", "Summits"])
        tbl.startBody()
        for contr in SC.contrasts:
            tbl.addRow([ Table.C(contr['test']),
                         Table.C(contr['control']),
                         Table.A(contr['macsPileup']),
                         Table.A(contr['macsNarrow']),
                         Table.A(contr['macsSummits']) ])
        tbl.toHTML(ACT.out)
        return True

### Homer

class HomerTags(Line):
    """Invoke Homer makeTagDirectory on BAM files for all conditions."""
    name = "Homer (tags)"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nhomer = 0
        for c in SC.conditions:
            cname = c['name']
            c['tags'] = cname + ".tags.d/"
            c['itags'] = cname + ".itags.d/"
            bams = SC.conditionBAMs(cname)
            LOG.log("BAMs for condition {}: {}", cname, bams)
            if not self.dry:
                ACT.submit("homer.qsub makeTagDirectory {} {}".format(c['tags'], " ".join(bams)), done="homer.@.done")
                nhomer += 1
            ibams = SC.conditionBAMs(cname, role='input')
            LOG.log("BAMs for condition {} (input): {}", cname, ibams)
            if not self.dry:
                ACT.submit("homer.qsub makeTagDirectory {} {}".format(c['itags'], " ".join(ibams)), done="homer.@.done")
                nhomer += 1
        return ACT.wait(("homer.@.done", nhomer))

class HomerPeaks(Line):
    """Invoke Homer findPeaks on a tag directory."""
    name = "Homer (peaks)"
    
    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        self.cnames = []
        self.peaks = []
        self.histones = []
        self.supers = []

        nhomer = 0
        for c in SC.conditions:
            self.cnames.append(c['name'])
            self.peaks.append(c['tags'] + "/peaks.txt")
            self.histones.append(c['tags'] + "/regions.txt")
            self.supers.append(c['tags'] + "/superEnhancers.txt")

            LOG.log("Finding peaks in {} (input: {})".format(c['tags'], c['itags']))
            if not self.dry:
                for mode in ["factor", "histone", "super"]:
                    ACT.submit("homer.qsub findPeaks {} {} {}".format(mode, c['tags'], c['itags']), done="homerp.@.done")
                    nhomer += 1
        ACT.wait(("homerp.@.done", nhomer))
        clscript = writeClassifyScript(SC.conditions, ACT.genesdb)
        if not self.dry:
            ACT.submit(clscript, done="homerc.@.done")
            ACT.wait(("homerc.@.done", 1))
        nhomer = 0
        for c in SC.conditions:
            if not self.dry:
                j = ACT.submit("homer.qsub annotate {} {}".format(c['tags'], "hg38"), done="homera.@.done")
                nhomer += 1
        return ACT.wait(("homera.@.done", nhomer))

    def PostExecute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        # Convert all peaks, regions, and enhancers files to Excel
        n = len(self.cnames)
        if ACT.missingOrStale("peaks.xlsx", other=self.peaks):
            args = ["{} -name {}".format(self.peaks[i], self.cnames[i]) for i in range(n) ]
            ACT.shell("module load dibig_tools; csvtoxls.py peaks.xlsx -q {}".format(" ".join(args)))
        if ACT.missingOrStale("regions.xlsx", other=self.histones):
            args = ["{} -name {}".format(self.histones[i], self.cnames[i]) for i in range(n) ]
            ACT.shell("module load dibig_tools; csvtoxls.py regions.xlsx -q {}".format(" ".join(args)))
        if ACT.missingOrStale("enhancers.xlsx", other=self.supers):
            args = ["{} -name {}".format(self.supers[i], self.cnames[i]) for i in range(n) ]
            ACT.shell("module load dibig_tools; csvtoxls.py enhancers.xlsx -q {}".format(" ".join(args)))

        # Generate plot of sites distribution
        classificationBarChart("peaks-classification.png", SC.conditions, "peaks", "Distribution of peak locations.")
        classificationBarChart("regions-classification.png", SC.conditions, "regions", "Distribution of region locations.")
        classificationBarChart("enhancers-classification.png", SC.conditions, "superEnhancers", "Distribution of enhancer locations.")
        return True

    def Report(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        ACT.scene("Peak finding")
        ACT.reportf("""Peak finding was performed using <b>HOMER</b>. Three detection types were performed: <i>Peaks</i>, <i>Regions</i>, 
and <i>SuperEnhancers</i>. Please see the <A href="http://homer.ucsd.edu/homer/ngs/peaks.html">HOMER</A> documentation for details. 
The following table provides links to the bedGraph files for each detection type in each condition, and to an Excel file containing 
the results of each detection type split by genomic region (e.g. Promoters, Exons, Introns, Intergenic). The <i>Combined</i> file
contains the results of each detection type in all conditions (one sheet per condition).""")
        tbl = Table.ScrollingTable(id='tblhomer', align='RL', caption="Results of peak detection with HOMER.")
        tbl.startHead()
        tbl.addHeaderRow(["Condition", "Count", "bedGraph", "Classified"])
        tbl.startBody()

        tbl.addSectionRow("Peaks")
        npeaks = 0
        for c in SC.conditions:
            path = c['tags'] + "peaks.bedGraph"
            nl = int(ACT.fileLines(path)) - 1
            npeaks += nl
            tbl.addRowHeader(c['name'])
            tbl.addRow([ Table.D(nl), Table.A(path), Table.A(c['peaksdistr']) ])
        tbl.addRowHeader("Combined")
        tbl.addRow( [ Table.D(npeaks), Table.A("peaks.xlsx"), "&nbsp;" ] )

        tbl.addSectionRow("Regions")
        npeaks = 0
        for c in SC.conditions:
            path = c['tags'] + "regions.bedGraph"
            nl = int(ACT.fileLines(path)) - 1
            npeaks += nl
            tbl.addRowHeader(c['name'])
            tbl.addRow([ Table.D(nl), Table.A(path), Table.A(c['regionsdistr']) ])
        tbl.addRowHeader("Combined")
        tbl.addRow( [ Table.D(npeaks), Table.A("regions.xlsx"), "&nbsp;" ] )

        tbl.addSectionRow("Enhancers")
        npeaks = 0
        for c in SC.conditions:
            path = c['tags'] + "superEnhancers.bedGraph"
            nl = int(ACT.fileLines(path)) - 1
            npeaks += nl
            tbl.addRowHeader(c['name'])
            tbl.addRow([ Table.D(nl), Table.A(path), Table.A(c['enhancersdistr']) ])
        tbl.addRowHeader("Combined")
        tbl.addRow( [ Table.D(npeaks), Table.A("enhancers.xlsx"), "&nbsp;" ] )

        tbl.toHTML(ACT.out)

        ACT.reportf("The following plots show the distribution of peak locations in the different conditions.")
        ACT.img("peaks-classification.png", thumbsize=200, description="Distribution of peak locations")
        ACT.img("regions-classification.png", thumbsize=200, description="Distribution of region locations")
        ACT.img("enhancers-classification.png", thumbsize=200, description="Distribution of enhancer locations")

        return True

class HomerDiffPeaks(Line):
    """Determine differential peaks."""
    name = "Homer (diff)"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nhomer = 0
        for contrast in SC.contrasts:
            ctrl = contrast['control']
            test = contrast['test']
            ctrlcond = SC.findCondition(ctrl)
            testcond = SC.findCondition(test)
            name1 = "{}.vs.{}.diffPeaks".format(test, ctrl)
            name2 = "{}.vs.{}.diffPeaks".format(ctrl, test)
            contrast['diff1'] = name1 + ".csv"
            contrast['diff2'] = name2 + ".csv"
            contrast['diffbedg1'] = name1 + ".bedGraph"
            contrast['diffbedg2'] = name2 + ".bedGraph"
            if not self.dry:
                ACT.submit("homer.qsub diffPeaks {} {} {}".format(testcond['tags'], ctrlcond['tags'], name1), done="homerd.@.done")
                ACT.submit("homer.qsub diffPeaks {} {} {}".format(ctrlcond['tags'], testcond['tags'], name2), done="homerd.@.done")
                nhomer += 2
        return ACT.wait(("homerd.@.done", nhomer))

    def PostExecute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc
        
        for contrast in SC.contrasts:
            ctrl = contrast['control']
            test = contrast['test']
            contrast['diffx'] = "{}.vs.{}.diffPeaks.xlsx".format(test, ctrl)
            if ACT.missingOrStale(contrast['diffx'], other=contrast['diff1']):
                ACT.shell("module load dibig_tools; csvtoxls.py {} -q {} {}".format(contrast['diffx'], contrast['diff1'], contrast['diff2']))
        return True

    def Report(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        ACT.scene("Differential peak finding")
        ACT.reportf("""Differential peak finding was performed using the HOMER <B>getDifferentialPeaks</B> command. The following files contain 
the peaks showing significant differences in each contrast.""")
        tbl = Table.ScrollingTable(id='tblhomerdiff', align='CCL', caption="Results of differential peak detection with HOMER.")
        tbl.startHead()
        tbl.addHeaderRow(["Test", "Control", "DiffPeaks"])
        tbl.startBody()
        for contr in SC.contrasts:
            tbl.addRow([ contr['test'],
                         contr['control'],
                         Table.A(contr['diffx']) ])
        tbl.toHTML(ACT.out)
        return True

class HomerMotifs(Line):
    """Perform motif identification on peaks."""
    name = "Homer (motifs)"
    genome = "hg38"

    def Setup(self):
        self.genome = Utils.dget('genome', self.properties, self.genome)
        return True

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        nhomer = 0
        for c in SC.conditions:
            outdir = c['name'] + ".motifs"
            peaks  = c['tags'] + "/peaks.txt"
            LOG.log("Finding motifs in {}, output to {}".format(peaks, outdir))
            # Add motfis dir to Zip file
            ACT._addToInclude(outdir + "/*")
            #ACT.shell("echo '{}/{}/*' >> .files".format(ACT.Name, outdir))

            if not self.dry:
                ACT.submit("homer.qsub findMotifs {} {} {}".format(peaks, self.genome, outdir), done="homerm.@.done")
                nhomer += 1
        return ACT.wait(("homerm.@.done", nhomer))

    def Report(self):
        ACT = self.actor
        SC = ACT.sc

        ACT.scene("Motif finding")
        ACT.reportf("""ChIP-Seq peaks were analyzed with the HOMER <b>findMotifs</b> function. The links in the following table lead to the motif finding report for each condition.""")
        tbl = Table.ScrollingTable(id='tblhomermotifs', align='CLL', caption="Results of motif finding with HOMER.")
        tbl.startHead()
        tbl.addHeaderRow(["Condition", "Known", "<i>de novo</i>"])
        tbl.startBody()
        for c in SC.conditions:
            outdir = c['name'] + ".motifs/"
            tbl.addRow([ c['name'], 
                         Table.A(outdir + "knownResults.html", target='homer'), 
                         Table.A(outdir + "homerResults.html", target='homer')])
        tbl.toHTML(ACT.out)
        return True

### Insert size analysis

class InsertSize(Line):
    """Compute the insert size distribution on the contents of a BAM file using Picard."""
    name = "Insert Size"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        npicard = 0
        for cond in SC.conditions:
            cond['insertsizecsv'] = cond['name'] + ".insert-size-metrics.csv"
            cond['insertsizepdf'] = cond['name'] + ".insert-size-metrics.pdf"
            cmdline = "picard.qsub CollectInsertSizeMetrics I={} O={} H={} M=0.5".format(cond['bam'], cond['insertsizecsv'], cond['insertsizepdf'])
            LOG.log("Executing: {}", cmdline)
            if not self.dry:
                ACT.submit(cmdline, done="isize.@.done")
                npicard += 1
        return ACT.wait(("isize.@.done", npicard))

    def Report(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        ACT.scene("Insert size distribution")
        ACT.reportf("""Insert size distribution was analyzed using the Picard <b>CollectInsertSizeMetrics</b> tool. The following PDF files contain the plot of insert size distribution in each condition.""")
        for cond in SC.conditions:
            ACT.file(cond['insertsizepdf'], description="Insert size distribution in {}".format(cond['name']))
        return True

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
                          Table.D(prereadcounts[sn]),
                          Table.D(readcounts[sn]), 
                          Table.D(prereadcounts[sn] - readcounts[sn]),
                          Table.P(1.0 * readcounts[sn] / prereadcounts[sn])])
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
        
        paired = SC.samples[0]['readsets'][0]['paired'] # Maybe this should be stored in SC?
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
            al = Utils.dget(sn, nalign, default=0) # number of aligned reads
            if paired:
                al = al/2
            tbl.addRow([ "<b>" + sn + "</b>",
                         Table.D(readcounts[sn]),
                         Table.D(al),
                         Table.P(al, Utils.dget(sn, readcounts, default=0)) ])
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

        nsplit = 0
        for smp in SC.samples:
            smp['bamC'] = smp['name'] + ".C.bam"
            smp['bamG'] = smp['name'] + ".G.bam"
            if ACT.missingOrStale(smp['bamC'], smp['bam']) or ACT.missingOrStale(smp['bamG'], smp['bam']):
                LOG.log("Splitting converted strands from {} into {} and {}.", smp['bam'], smp['bamC'], smp['bamG'])
                if not self.dry:
                    ACT.submit("generic.qsub bisconv.py split {} {} {} module:dibig_tools".format(smp['bam'], smp['bamC'], smp['bamG']), done="split.@.done")
                    nsplit += 1
        ACT.wait(("split.@.done", nsplit))

        ncall = 0
        for cond in SC.conditions:
            cbams = ",".join([ smp['bamC'] for smp in SC.conditionSamples(cond['name']) ])
            gbams = ",".join([ smp['bamG'] for smp in SC.conditionSamples(cond['name']) ])
            cond['callbed'] = cond['name'] + ".bed"
            cond['callbedg'] = cond['name'] + ".bedGraph"
            cond['callreport'] = cond['name'] + "-report.csv"
            cond['lcreport'] = cond['name'] + "-lc.csv"
            cond['histreport'] = cond['name'] + "-hist.csv"
            cond['matreport'] = cond['name'] + "-mat.csv"
            cond['winavg'] = cond['name'] + "-winavg.csv"
            cond['winmat'] = cond['name'] + "-winmat.csv"
            cond['winavgx'] = cond['name'] + "-winavg.xlsx"
            cond['winmatx'] = cond['name'] + "-winmat.xlsx"

            if ACT.strand == None:
                strand = ""
            else:
                strand = "strand=" + ACT.strand
            LOG.log("Calling methylation on BAM files for sample {}.".format(cond['name']))
            if not self.dry:
                ACT.submit("cscall.qsub {} {} {} {} out={} report={} lc={} hist={} mat={} mind={} mins={} {}".format(cbams, gbams, ACT.siteindex, ACT.reference, cond['callbed'], cond['callreport'], cond['lcreport'], cond['histreport'], cond['matreport'], ACT.mindepth, ACT.minsamples, strand), done="mcall.@.done")
                ncall += 1
        return ACT.wait(("mcall.@.done", ncall))

    def PostExecute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        allwinavg = []
        for cond in SC.conditions:
            if ACT.missingOrStale(cond['callbedg'], cond['callbed']):
                ACT.shell('echo "track type=bedGraph" > {}; cut -f 1,2,3,4 {} >> {}'.format(cond['callbedg'], cond['callbed'], cond['callbedg']))
            allwinavg.append(cond['winavg'])
            if ACT.missingOrStale(cond['winavg'], cond['callbed']):
                ACT.shell("module load dibig_tools; dmaptools.py winavg -o {} {}", cond['winavg'], cond['callbed'])
            if ACT.missingOrStale(cond['winavgx'], cond['winavg']):
                ACT.shell("module load dibig_tools; csvtoxls.py {} -q {}", cond['winavgx'], cond['winavg'])
            if ACT.missingOrStale(cond['winmat'], cond['matreport']):
                ACT.shell("module load dibig_tools; dmaptools.py winmat -o {} {}", cond['winmat'], cond['matreport'])
            if ACT.missingOrStale(cond['winmatx'], cond['winmat']):
                ACT.shell("module load dibig_tools; csvtoxls.py {} -q {}", cond['winmatx'], cond['winmat'])
            cond['callreportx'] = cond['name'] + "-report.xlsx"
            ACT.shell("module load dibig_tools; csvtoxls.py {} {} -firstrowhdr -name Report {} -firstrowhdr -name Histogram {} -firstrowhdr -name LoneC",
                      cond['callreportx'], cond['callreport'], cond['histreport'], cond['lcreport'])

        # Combine all winavg files together
        if ACT.missingOrStale("Merged.winavg.csv", allwinavg):
            ACT.shell("module load dibig_tools; dmaptools.py cmerge -o Merged.winavg.csv {}; csvtoxls.py Merged.winavg.xlsx -q Merged.winavg.csv".format(" ".join(allwinavg)))

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
        ACT.reportf("""The following table provides a summary of the number of sites identified for each experimental condition. Follow the links to the reports for more detailed information.
The bedGraph file contains the % methylation level at each detected site, and is suitable for IGV or the UCSC Genome Browser.""")
        tbx = Table.ScrollingTable(id='tblcscall', align="LRRRRCCCC", caption="Summary of methylation calling.")
        tbx.startHead()
        tbx.addHeaderRow(["Condition", "Sites", "Total basepairs", "Total site coverage", "Avg site coverage", "Full report", "bedGraph"])
        tbx.startBody()
        for cond in SC.conditions:
            data = getCscallReport(cond)
            if data:
                tbx.addRow([ "<b>" + cond['name'] + "</b>",
                             Table.D(int(data[1])),
                             Table.D(int(data[0])),
                             Table.D(int(data[2])),
                             Table.F(data[3]),
                             Table.A(cond['callreportx']),
                             Table.A(cond['callbedg']) ])
        tbx.toHTML(ACT.out)

        ACT.reportf("""The following table provides links to additional reports on methylation calling. The WinAvg files contain the average
methylation in each window (averaged over all sites contained in the window and all replicates). The Raw Data files contain the % methylation
level at each detected site in each replicate, preceded by their average and standard deviation.""")
        tbx = Table.ScrollingTable(id='tblcscall2', align="LRRRRCCCC", caption="Additional methylation reports.")
        tbx.startHead()
        tbx.addHeaderRow(["Condition", "WinAvg (CSV)", "WinAvg (Excel)", "Raw data (CSV)", "Plots"])
        tbx.startBody()
        for cond in SC.conditions:
            if data:
                tbx.addRow([ Table.C("<b>" + cond['name'] + "</b>"),
                             Table.A(cond['winavg']),
                             Table.A(cond['winavgx']),
                             Table.A(cond['matreport']),
                             Table.A(cond['methplots']) ])
        tbx.toHTML(ACT.out)

        ACT.file("Merged.winavg.xlsx", description="Combined winavg file.")
        return True
        
### MCOMP

class MCOMP(Line):
    name = "MCOMP"

    def Execute(self):
        ACT = self.actor
        SC  = ACT.sc
        LOG = ACT.log

        if len(SC.contrasts) > 0:
            ncomp = 0
            subs = writeMCOMPqsub(ACT)
            if not self.dry:
                for sub in subs:
                    ncomp += 1
                    ACT.submit(sub, done="mcomp.@.done")
            ACT.wait(("mcomp.@.done", ncomp))
            LOG.log("Computing differential methylation summary.")
            self.prefix = Utils.id_generator(10, prefix='tmp-')
            ACT.diffmeth = diffmethSummary(ACT, self.prefix, self.dry)
        return True

    def Report(self):
        ACT = self.actor
        SC  = ACT.sc
        LOG = ACT.log

        if len(SC.contrasts) == 0:
            return True

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
                         Table.D(d['total']),
                         "{} ({})".format(totsig, Utils.pct(totsig, d['total'])),
                         "{} ({})".format(Utils.UP(Table.D(d['++'])), Utils.pct(d['++'], totsig)),
                         "{} ({})".format(Utils.up(Table.D(d['+'])), Utils.pct(d['+'], totsig)),
                         "{} ({})".format(Utils.down(Table.D(d['-'])), Utils.pct(d['-'], totsig)),
                         "{} ({})".format(Utils.DOWN(Table.D(d['--'])), Utils.pct(d['--'], totsig))
                         ])

        tbl.toHTML(ACT.out)

        ACT.reportf("""<BR><BR>The following table reports the <b>average genome-wide metyhlation</b> in all replicates of each condition. 
The P-value column reports the P-value of the significance of the difference between the two sets of values, computed using a 
<b>two-sample t-test</b>. The last column indicates whether the P-value is below <b>1%</b>.""")
        tbl = Table.ScrollingTable(id='tblmcomp3', align='CCRRRC', caption="Difference in average genome-wide methylation between samples.")
        tbl.startHead()
        tbl.addHeaderRow(["Test", "Control", "Test Avgs", "Control Avgs", "P-value", "Sig"])
        tbl.startBody()
        for contr in SC.contrasts:
            data = runAvgMeth(ACT, contr['test'], contr['control'])
            tbl.addRow(["<b>" + contr['test'] + "</b>",
                        "<b>" + contr['control'] + "</b>",
                        "<BR>".join(data["MethRates1"].split(",")),
                        "<BR>".join(data["MethRates2"].split(",")),
                        data["P-value"],
                        data["Significant"]])
        tbl.toHTML(ACT.out)

        ACT.reportf("""<br><br>The following table describes the difference in the distribution of methylated sites between conditions. The
results are also available as Excel file {}.""", PipelineUtils.linkify("methylHist.xlsx", "methylHist.xlsx"))
        outfiles = []
        tbl = Table.ScrollingTable(id='tblmcomp4', align='CRRRC', caption="Difference in distribution of methylation levels between samples.")
        tbl.startHead()
        tbl.addHeaderRow(["Bin", "Test Pct", "Control Pct", "P-value", "Sig"])
        tbl.startBody()
        for contr in SC.contrasts:
            test = contr['test']
            ctrl = contr['control']
            label = test + " vs " + ctrl
            outfile = test + ".vs." + ctrl + ".csv"
            outfiles.append(outfile)
            tbl.addSectionRow(label)
            result = runMethHist(ACT, test, ctrl, outfile)
            for row in result:
                tbl.addRow([row[0], Utils.pct(float(row[3]), 1), Utils.pct(float(row[4]), 1), row[1], row[2]])
        tbl.toHTML(ACT.out)
        ACT.shell("module load dibig_tools; csvtoxls.py -q methylHist.xlsx " + " ".join(outfiles))

        ACT.reportf("""<br><br>The following table provides links to files containing detailed differential methylation data. The <i>Significant</i> files contain sites showing differential methylation only, while the <i>Full</i> files contain data for all sites. The <i>BED</i> files can be loaded into a genome browser (e.g., IGV), and the <i>Plots</i> files contain plots of differential methylation across chromosomes.""")

        tbl = Table.ScrollingTable(id='tblmcomp2', align="CCCCCC", caption="Tables of differential methylation data.")
        tbl.startHead()
        tbl.addHeaderRow(["Test", "Control", "Significant", "Full", "BED", "Plots"])
        tbl.startBody()
        for d in ACT.diffmeth:
            tbl.addRow([ "<b>" + d['test'] + "</b>",
                         "<b>" + d['control'] + "</b>",
                         Table.A(d['sigfilex']),
                         Table.A(d['fullfilex']),
                         Table.A(d['bedfile']),
                         Table.A(d['dmcfile']) ])

        tbl.toHTML(ACT.out)
        return True

### DMR

class DMR(Line):
    name = "Differentially methylated regions"

    def makeCmdline(self, ACT, infile1, infile2, outfile):
        return "generic.qsub dmaptools.py dmr -o {} -w {} -s {} -c {} -d {} -p {} -g {} {} {} module:dibig_tools".format(
            outfile, ACT.dmr_winsize, ACT.dmr_minsites, ACT.dmr_mincov, ACT.dmr_diff, ACT.dmr_pval, ACT.dmr_gap, infile1, infile2)

    def Execute(self):
        ACT = self.actor
        SC  = ACT.sc
        LOG = ACT.log
        
        ndmr = 0
        for contr in SC.contrasts:
            test = contr['test']
            ctrl = contr['control']
            testcond = SC.findCondition(test)
            ctrlcond = SC.findCondition(ctrl)
            outfile = "{}.vs.{}.dmr".format(test, ctrl)
            contr['dmrfile'] = outfile + ".csv"
            contr['dmrxls'] = outfile + ".xlsx"
            cmdline = self.makeCmdline(ACT, testcond['callbed'], ctrlcond['callbed'], contr['dmrfile'])
            LOG.log("Detecting DMRs in {} vs {}", test, ctrl)
            LOG.log("Commandline: {}", cmdline)
            if not self.dry:
                ACT.submit(cmdline, done="dmr.@.done")
                ndmr += 1
        return ACT.wait(("dmr.@.done", ndmr))

    def PostExecute(self):
        ACT = self.actor
        SC  = ACT.sc

        alldmrs = []
        for contr in SC.contrasts:
            alldmrs.append(contr['dmrfile'])
            ACT.shell("module load dibig_tools; csvtoxls.py -q {} {}".format(contr['dmrxls'], contr['dmrfile']))

        if ACT.missingOrStale("Merged.dmr.csv", alldmrs):
            ACT.shell("module load dibig_tools; dmaptools.py cmerge -o Merged.dmr.csv {}; csvtoxls.py Merged.dmr.xlsx -q Merged.dmr.csv".format(" ".join(alldmrs)))
        return True

    def Report(self):
        ACT = self.actor
        SC  = ACT.sc

        ACT.scene("Differentially methylated regions")
        ACT.reportf("""The following table reports the differentially methylated regions (DMRs) detected in each contrast, 
using the method described in {}.""".format(PipelineUtils.linkify("http://www.sciencedirect.com/science/article/pii/S0092867412014304", "Stroud et al., Cell 2013")))
        tbl = Table.ScrollingTable(id='tbldmr', align="CCRC", caption="Tables of differentially methylated regions (DMRs).")
        tbl.startHead()
        tbl.addHeaderRow(["Test", "Control", "DMRs", "DMR file (CSV)", "DMR file (Excel)"])
        tbl.startBody()
        for contr in SC.contrasts:
            tbl.addRow(["<b>" + contr['test'] + "</b>",
                        "<b>" + contr['control'] + "</b>",
                        int(ACT.fileLines(contr['dmrfile'])) - 1,
                        Table.A(contr['dmrfile']),
                        Table.A(contr['dmrxls'])])
        tbl.toHTML(ACT.out)

        ACT.file("Merged.dmr.xlsx", "Merged DMR file.")
        return True

### GeneMeth

class GeneMeth(Line):
    name = "Gene methylation"

    def makeCmdline(self, ACT, infile, outfile, idx=0, srt=False):
        cmdline = "generic.qsub genediffmeth "
        if ACT.regions:
            cmdline += " -r " + ACT.regions[idx]
        if ACT.size:
            cmdline += " -d " + ACT.size[idx]
        if ACT.mode:
            cmdline += " -m " + ACT.mode[idx]
        if ACT.nsites:
            cmdline += " -l " + ACT.nsites[idx]
        if srt:
            cmdline += " -s "
        return cmdline + " " + infile + " " + ACT.genesdb + " " + outfile + " module:dibig_tools"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        
        ngd = 0   # Number of genediffmeth jobs submitted
        l = 0 # Counter for different genediffmeth runs (with different parameter sets)
        runlabel = ""

        if ACT.genesdb:
            for l in range(len(ACT.regions)):
                runlabel = "{:02d}".format(l+1)
                #print "l={}, runlabel={}".format(l, runlabel)
                #raw_input()
                for dm in ACT.diffmeth:
                    outfile = "{}.vs.{}.genes{}.sig.csv".format(dm['test'], dm['control'], runlabel)
                    dm['genesigfile'+runlabel] = outfile
                    cmdline = self.makeCmdline(ACT, dm['sigfile'], outfile, idx=l, srt=True)
                    LOG.log("Generating differential methylation file for genes containing significant sites.")
                    LOG.log("Commandline: {}", cmdline)
                    if not self.dry:
                        ACT.submit(cmdline, done="gdm.@.done")
                        ngd += 1
                    outfile = "{}.vs.{}.genes{}.full.csv".format(dm['test'], dm['control'], runlabel)
                    dm['genefullfile'+runlabel] = outfile
                    cmdline = self.makeCmdline(ACT, dm['fullfile'], outfile, idx=l)
                    LOG.log("Generating differential methylation file for all genes.")
                    LOG.log("Commandline: {}", cmdline)
                    if not self.dry:
                        ACT.submit(cmdline, done="gdm.@.done")
                        ngd += 1
        return ACT.wait(("gdm.@.done", ngd))

    def PostExecute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        l = 0 # Counter for different genediffmeth runs (with different parameter sets)
        rl = ""

        for dm in ACT.diffmeth:
            for l in range(len(ACT.regions)):
                rl = "{:02d}".format(l+1)
                dm['genesigfilex'+rl] = "{}.vs.{}.genes{}.sig.xlsx".format(dm['test'], dm['control'], rl)
                dm['genefullfilex'+rl] = "{}.vs.{}.genes{}.full.xlsx".format(dm['test'], dm['control'], rl)
                if ACT.missingOrStale(dm['genesigfilex'+rl], dm['genesigfile'+rl]):
                    ACT.shell("module load dibig_tools; csvtoxls.py {} -q {};",
                              dm['genesigfilex'+rl], dm['genesigfile'+rl])
                if ACT.missingOrStale(dm['genefullfilex'+rl], dm['genefullfile'+rl]):
                    ACT.shell("module load dibig_tools; csvtoxls.py {} -q {};",
                              dm['genefullfilex'+rl], dm['genefullfile'+rl])
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

    def printRegnames(self, r):
        regnames = {'p': 'promoter', 'b': 'gene body', 'e': 'exons', 'i': 'introns', 'd': 'downstream'}
        words = [ regnames[w] for w in r ]
        return ", ".join(words)

    def Report(self):
        ACT = self.actor
        LOG = ACT.log

        l = 0 # Counter for different genediffmeth runs (with different parameter sets)
        rl = ""

        if ACT.genesdb:
            ACT.scene("Gene methylation analysis")
            ACT.reportf("""<P>The following table provides links to files containing information about methylation in genes. For each contrast, two files are provided, one 
    containing only genes showing evidence of differential methylation (ordered from highest hypermethylation to highest hypomethylation) and one containing all
    genes (ordered by chromosome and position).</P>""")

            for l in range(len(ACT.regions)):
                rl = "{:02d}".format(l+1)
                ACT.reportf("""<P>Region: <b>{}</b>; Upstream/downstream: <b>{}</b>; Operator: <b>{}</b>; Number of sites: <b>{}</b>.</P>""".format(
                        self.printRegnames(ACT.regions[l]), ACT.size[l], ACT.mode[l], ACT.nsites[l]))
                
                tbl = Table.ScrollingTable(id='tblgmeth'+rl, align="CCRRRRRCC", caption="Differential methylation in genes. For each contrast, the table lists the number of genes showing evidence of differential methylation, and the breakdown into the four differential methylation classes.")
                tbl.startHead()
                tbl.addHeaderRow(["Test", "Control", "Genes", "++", "+", "-", "--", "Significant", "Full"])
                tbl.startBody()
                for d in ACT.diffmeth:
                    counts = self.countSigGenes(d['genesigfile'+rl], ACT.diff)
                    print counts
                    tot = counts['total']
                    tbl.addRow([ "<b>" + d['test'] + "</b>",
                                 "<b>" + d['control'] + "</b>",
                                 Table.D(tot),
                                 "{} ({})".format(Utils.UP(Table.D(counts['++'])),   Utils.pct(counts['++'], tot)),
                                 "{} ({})".format(Utils.up(Table.D(counts['+'])),    Utils.pct(counts['+'], tot)),
                                 "{} ({})".format(Utils.down(Table.D(counts['-'])),  Utils.pct(counts['-'], tot)),
                                 "{} ({})".format(Utils.DOWN(Table.D(counts['--'])), Utils.pct(counts['--'], tot)),
                                 Table.A(d['genesigfilex'+rl]),
                                 Table.A(d['genefullfilex'+rl]) ])

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

### SNP calling

class SNPcallGATK(Line):
    name = "SNP calling with GATK"

    def Verify(self):
        """Verify that the reference is indexed and has a dictionary."""
        ACT = self.actor
        
        ref = ACT.reference
        refidx = ref + ".fai"
        refdict = ACT.setFileExt(ref, ".dict")
        if ACT.missingOrStale(refidx, other=ref):
            return "The reference file needs to be indexed with samtools ({} not found our out of date).".format(refidx)
        if ACT.missingOrStale(refdict, other=ref):
            return "The reference file needs a Picard dictionary ({} not found our out of date).".format(refdict)
        return True

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC  = ACT.sc

        ngatk = 0
        LOG.log
        for smp in SC.samples:
            smp['gvcf'] = smp['name'] + ".g.vcf"
            ACT.exclude(smp['gvcf'])
            if not self.dry:
                ACT.submit("gatk-snpcall-gvcf.qsub {} {} {}".format(ACT.reference, smp['bam'], smp['gvcf']), done="gvcf.@.done")
                ngatk += 1
                
        ACT.wait(("gvcf.@.done", ngatk))

        ngatk = 0
        for cond in SC.conditions:
            cond['vcf'] = cond['name'] + ".vcf"
            tmpfile = "tmp-" + cond['name'] + "-gvcf.txt"
            with open(tmpfile, "w") as out:
                for smp in SC.conditionSamples(cond):
                    out.write(smp['gvcf'] + "\n")
            if not self.dry:
                ACT.submit("gvcf.qsub {} {} {}".format(ACT.reference, tmpfile, cond['vcf']), done="gvcf.@.done")
                ngatk += 1
        ACT.wait(("gvcf.@.done", ngatk))

        # Clean the VCF files (get rid of bad chroms).
        for cond in SC.conditions:
            cond['origvcf'] = cond['name'] + ".pre.vcf"
            if ACT.missingOrStale(cond['origvcf'], other=cond['vcf']):
                ACT.copy(cond['vcf'], cond['origvcf'])
            ACT.shell("grep -v chrUn {} | grep -v _random > {}".format(cond['origvcf'], cond['vcf']))

        return True

    def PostExecute(self):
        ACT = self.actor
        LOG = ACT.log
        SC  = ACT.sc

        # Run snpEff on VCF for each condition
        nsnpeff = 0
        if ACT.snpEffDB:
            for cond in SC.conditions:
                cond['snpeff'] = snpeffRec(cond['vcf'])
                if not self.dry:
                    cond['snpeff'].runsnpEff(ACT)
                    nsnpeff += 1
        ACT.wait(("snpeff.@.done", nsnpeff))

        # Convert snpEff reports
        for cond in SC.conditions:
            cond['snpeff'].countSNPs(ACT)
            cond['snpeff'].toExcel(ACT)
            cond['snpeff'].getFiltered(ACT, xls=cond['name'] + ".filtered.xlsx")

        return True

    def Report(self):
        ACT = self.actor
        LOG = ACT.log
        SC  = ACT.sc

        ACT.scene("Variant calling")
        ACT.reportf("""Variant calling was performed using the GATK pipeline. Variants were then annotated using snpEff.
The following table provides links to the main VCF output file, the snpEff report and a table showing the number of variants
in each functional class for each transcript.""")

        tbl = Table.ScrollingTable(id='tblsnps', caption="VCF and annotation files for each condition.")
        tbl.startHead()
        tbl.addHeaderRow(["Condition", "Variants", "VCF", "High impact", "snpEff report", "snpEff genes"])
        tbl.startBody()
        for cond in SC.conditions:
            tbl.addRow([ Table.C(cond['name']),
                         Table.D(cond['snpeff'].nsnps),
                         Table.A(cond['vcf']),
                         Table.A(cond['snpeff'].high),
                         Table.A(cond['snpeff'].html),
                         Table.A(cond['snpeff'].xlsx) ])
        tbl.toHTML(ACT.out)

        ACT.reportf("""The following table reports the number of SNPs in the specified functional classes. Full details about the SNPs can be found in the Excel files below the table.""")

        tbl = Table.ScrollingTable(id='tblsnpannot', caption="Subsets of SNPs in various functional classes.")
        tbl.startHead()
        tbl.addHeaderRow(["Condition", "Effect", "Variants", "VCF"])
        tbl.startBody()
        for cond in SC.conditions:
            filt = cond['snpeff'].filtered
            nf = sum([1 if f[4] > 0 else 0 for f in filt])
            tbl.addRowHeader(cond['name'], rowspan=nf)
            for f in filt:
                if f[4] > 0:
                    tbl.addRow([ Table.C(f[0]), Table.D(f[4]), Table.A(f[1]) ])
        tbl.toHTML(ACT.out)

        for cond in SC.conditions:
            ACT.file(cond['name'] + ".filtered.xlsx", description="SNPs with functional annotations in {}.".format(cond['name']))

        return True

class SNPdiff(Line):
    """For each contrast, determine (a) common SNPs, (b) SNPs only in first condition, (c) SNPs only in second conditions."""

    name = "Common and unique SNPs in contrasts."

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC  = ACT.sc

        # Now do the diff
        for contr in SC.contrasts:
            test = contr['test']
            ctrl = contr['control']
            testcond = SC.findCondition(test)
            ctrlcond = SC.findCondition(ctrl)
            out = "{}.vs.{}".format(ctrl, test)
            contr['sitediff'] = out + ".diff.sites_in_files"
            ACT.shell("module load vcftools; vcftools --vcf {} --diff {} --diff-site --out {}".format(ctrlcond['vcf'], testcond['vcf'], out))
            # Now we have one 'sitediff' file for each contrast. Split the VCF files based on it.
            contr['splitvcf'] = splitVCF(out, contr['sitediff'], ctrlcond['vcf'], testcond['vcf'])
        return True

    def PostExecute(self):
        ACT = self.actor
        LOG = ACT.log
        SC  = ACT.sc

        # Run snpEff on the three VCFs for each condition
        if ACT.snpEffDB:
            nsnpeff = 0
            for contr in SC.contrasts:
                contr['snpeff'] = [ snpeffRec(v) for v in contr['splitvcf'] ]
                if not self.dry:
                    for sr in contr['snpeff']:
                        sr.runsnpEff(ACT)
                        nsnpeff += 1

            ACT.wait(("snpeff.@.done", nsnpeff))

            # Convert snpEff reports
            for contr in SC.contrasts:
                for sr, tag in zip(contr['snpeff'], ['A', 'B', 'both']):
                    sr.countSNPs(ACT)
                    sr.toExcel(ACT)
                    sr.getFiltered(ACT, xls="{}.vs.{}.{}.filtered.xlsx".format(contr['test'], contr['control'], tag))
        return True

    def Report(self):
        ACT = self.actor
        LOG = ACT.log
        SC  = ACT.sc

        ACT.scene("Differential variant calling")
        ACT.reportf("""For each contrast, the following table reports the SNPs present in the first file only, second file only, or both files.""")

        tbl = Table.ScrollingTable(id='tbldsnps', caption="VCF and annotation files for each condition.")
        tbl.startHead()
        tbl.addHeaderRow(["Contrast", "SNP set", "Variants", "VCF", "snpEff report", "snpEff genes"])
        tbl.startBody()
        for contr in SC.contrasts:
            [vcf1, vcf2, vcfboth] = contr['snpeff']
            test = contr['test']
            ctrl = contr['control']
            tbl.addRowHeader("{} vs {}".format(test, ctrl), rowspan=3)
            tbl.addRow([ "Ctrl only", Table.D(vcf1.nsnps), Table.A(vcf1.snpeff), Table.A(vcf1.html), Table.A(vcf1.xlsx)])
            tbl.addRow([ "Common",    Table.D(vcfboth.nsnps), Table.A(vcfboth.snpeff), Table.A(vcfboth.html), Table.A(vcfboth.xlsx)])
            tbl.addRow([ "Test only", Table.D(vcf2.nsnps), Table.A(vcf2.snpeff), Table.A(vcf2.html), Table.A(vcf2.xlsx)])

        tbl.toHTML(ACT.out)

        ACT.reportf("""The following table reports the number of SNPs in the specified functional classes. Full details about the SNPs can be found in the Excel files below the table.""")

        tbl = Table.ScrollingTable(id='tbldsnpannot', caption="Subsets of SNPs in various functional classes.")
        tbl.startHead()
        tbl.addHeaderRow(["Condition", "Effect", "Variants", "VCF"])
        tbl.startBody()
        for contr in SC.contrasts:
            [vcf1, vcf2, vcfboth] = contr['snpeff']
            test = contr['test']
            ctrl = contr['control']

            filt = vcf1.filtered
            nf = sum([1 if f[4] > 0 else 0 for f in filt])
            tbl.addRowHeader("{} vs {}<BR>Ctrl only".format(test, ctrl), rowspan=nf)
            for f in filt:
                if f[4] > 0:
                    tbl.addRow([ Table.C(f[0]), Table.D(f[4]), Table.A(f[1]) ])

            filt = vcfboth.filtered
            nf = sum([1 if f[4] > 0 else 0 for f in filt])
            tbl.addRowHeader("{} vs {}<BR>Common".format(test, ctrl), rowspan=nf)
            for f in filt:
                if f[4] > 0:
                    tbl.addRow([ Table.C(f[0]), Table.D(f[4]), Table.A(f[1]) ])

            filt = vcf2.filtered
            nf = sum([1 if f[4] > 0 else 0 for f in filt])
            tbl.addRowHeader("{} vs {}<BR>Test only".format(test, ctrl), rowspan=nf)
            for f in filt:
                if f[4] > 0:
                    tbl.addRow([ Table.C(f[0]), Table.D(f[4]), Table.A(f[1]) ])

        tbl.toHTML(ACT.out)

        for contr in SC.contrasts:
            test = contr['test']
            ctrl = contr['control']
            ACT.file("{}.vs.{}.{}.filtered.xlsx".format(contr['test'], contr['control'], 'A'))
            ACT.file("{}.vs.{}.{}.filtered.xlsx".format(contr['test'], contr['control'], 'B'))
            ACT.file("{}.vs.{}.{}.filtered.xlsx".format(contr['test'], contr['control'], 'both'))

        return True

### ATAC

class atacQC(Line):
    name = "ATAC quality control"

    # This should be moved to atacutils
    def convertHmap(self, filename, outfile):

        def getAvg(filename):
            tot = 0.0
            nrows = 0
            with open(filename, "r") as f:
                for line in f:
                    data = line.rstrip("\r\n").split("\t")[5:]
                    for i in range(100):
                        tot += float(data[i])
                    for i in range(1900, 2000):
                        tot += float(data[i])
                    nrows += 1
            return tot / nrows / 200

        base = getAvg(filename)
        # sys.stderr.write("Base: {}\n".format(base))
        vec = [0.0]*2000
        nrows = 0
        with open(filename, "r") as f:
            hdr = f.readline()
            for line in f:
                nrows += 1
                data = line.rstrip("\r\n").split("\t")[5:]
                for i in range(2000):
                    vec[i] += float(data[i])

        with open(outfile, "w") as out:
            for i in range(2000):
                out.write("{}\t{}\n".format(i - 1000, vec[i] / nrows / base))

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC  = ACT.sc

        ## Compute frip scores
        nfrip = 0
        for cond in SC.conditions:
            if not self.dry:
                ACT.submit("frip.qsub {} {} {}".format(cond['bam'], cond['macsdir'] + "/" + cond['name'] + "_peaks.csv", cond['name'] + ".frip.txt"), done="frip.@.done")
                nfrip += 1

        regs = "gene-regions.bed"
        if ACT.missingOrStale(regs):
            ACT.shell("module load dibig_tools; genes.py region -db {} -r u -d 1000 -f t all | grep -v chrMT > {}".format(ACT.genesdb, regs))

        nregs = 0
        for cond in SC.conditions:
            cond['hmap'] = cond['name'] + ".hmap.csv"
            cmdline = "generic.qsub regionsCount.py -w 2000 -o {} {} {} module:dibig_tools".format(cond['hmap'], cond['bam'], regs)
            if not self.dry:
                ACT.submit(cmdline, done="nregs.@.done")
                nregs += 1

        ACT.wait(("nregs.@.done", nregs), ("nfrip.@.done", nfrip))

        # Generate plots with gnuplot
        with open("enrich.gp", "w") as out:
            out.write("set term png size 1000,600\n")
            for cond in SC.conditions:
                cond['tssenrichpng'] = cond['name'] + ".tssenrich.png"
                cond['hmapcsv'] = cond['name'] + ".enscore.csv"
                if ACT.missingOrStale(cond['hmapcsv'], other=cond['hmap']):
                    self.convertHmap(cond['hmap'], cond['hmapcsv'])
                out.write("""
set output "{}"
set grid
set title "{} - TSS Enrichment Score"
set xlabel "Position"
set ylabel "Enrichment"
# set yrange [0.9:2]
plot "{}" with lines
""".format(cond['tssenrichpng'], cond['name'], cond['hmapcsv']))

        ACT.shell("module load gnuplot; gnuplot enrich.gp")
        return True

    def Report(self):
        ACT = self.actor
        LOG = ACT.log
        SC  = ACT.sc
        out = ACT.out

        ACT.scene("ATAC Quality Control")
        ACT.reportf("The following plots display TSS enrichment for each condition.")
        out.write("""<CENTER><TABLE class='figure'><TR><TD class='figure'>\n""")
        out.write("""  <TABLE width='100%'><TR><TD align='center' valign='top'>\n""")
        out.write("""    <TABLE class='legend'><TH class='legend'>Condition</TH>\n""")
        img1 = SC.conditions[0]['tssenrichpng']
        for cond in SC.conditions:
            out.write("""    <TD onclick="document.getElementById('tssenrich').src='{}';" class='legend'>{}</TD>\n""".format(cond['tssenrichpng'], cond['name']))
        out.write("""    </TR></TABLE>
  </TD></TR><TR><TD align='center' valign='top'>
    <A href='{}'><IMG class='bordered' id='tssenrich' src='{}' width='800'/></A>
  </TD></TR></TABLE>
</TD></TR>
</TABLE></CENTER>
""".format(img1, img1))

### UCSC hubs

class UCSCHub(Line):
    name = "Hub"
    kind = "rnaseq"
    hubname = ""
    dirname = "hub"

    def Setup(self):
        ACT = self.actor
        self.kind = Utils.dget('kind', self.properties, self.kind)
        self.hubname = ACT.getConf("hubname", "Hub") or ACT.title
        self.dirname = ACT.getConf("dirname", "Hub") or ACT.dirname

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc
        u = Utils.UCSCHub(self.hubname, ACT.getConf("shortLabel", "Hub"), ACT.getConf("longLabel", "Hub"), 
                          genome=ACT.getConf("genome", "Hub"), 
                          sizes=ACT.getConf("sizes", "Hub"), 
                          dirname=self.dirname, 
                          email=ACT.getConf("email", "Hub"))
        u.generate()
        if self.kind == 'rnaseq':
            u.startContainer("samples", "WIG tracks for samples", "WIGsamples", "bigWig")
            for smp in SC.samples:
                u.addWig(smp['wigfile'], smp['name'], smp['name'], "{} rnaseq", alwaysZero="on", smoothingWindow="4", visibility="dense")
            u.endContainer()
            u.startContainer("conditions", "WIG tracks for conditions", "WIGconditions", "bigWig")
            for cond in SC.conditions:
                u.addWig(cond['wigfile'], cond['name'], cond['name'], "{} rnaseq", alwaysZero="on", smoothingWindow="4", visibility="dense")
            u.endContainer()

        elif self.kind == 'chipseq':
            u.startContainer("samples", "WIG (by sample)", "WIGsamples", "bigWig")
            for smp in SC.samples:
                u.addWig(smp['wigfile'], smp['name'], smp['name'], "{} chipseq", alwaysZero="on", smoothingWindow="4", visibility="dense")
            u.endContainer()
            u.startContainer("conditions", "WIG (by condition)", "WIGconditions", "bigWig")
            for cond in SC.conditions:
                u.addWig(cond['wigfile'], cond['name'], cond['name'], "{} chipseq", alwaysZero="on", smoothingWindow="4", visibility="dense")
            u.endContainer()
            u.startContainer("pkbdg", "peaks", "BDGpeaks", "bigWig")
            for cond in SC.conditions:
                u.addBedGraph(cond['tags'] + "peaks.bedGraph", cond['name'] + "-peaks", cond['name'] + "-peaks", "{} - chipseq")
            u.endContainer()
            u.startContainer("rgbdg", "regions", "BDGregions", "bigWig")
            for cond in SC.conditions:
                u.addBedGraph(cond['tags'] + "regions.bedGraph", cond['name'] + "-regions", cond['name'] + "-regions", "{} - chipseq")
            u.endContainer()
            u.startContainer("sebdg", "enhancers", "BDGenhancers", "bigWig")
            for cond in SC.conditions:
                u.addBedGraph(cond['tags'] + "superEnhancers.bedGraph", cond['name'] + "-enhancers", cond['name'] + "-enhancers", "{} - chipseq")
            u.endContainer()
            u.startContainer("diffbdg", "differential peaks", "BDGdiffpeaks", "bigWig")
            for contr in SC.contrasts:
                u.addBedGraph(contr['diffbedg1'], contr['name'] + "-diffPeaks1", contr['name'] + "-diffPeaks1", "{} - chipseq")
                u.addBedGraph(contr['diffbedg2'], contr['name'] + "-diffPeaks2", contr['name'] + "-diffPeaks2", "{} - chipseq")
            u.endContainer

        elif self.kind == 'atacseq':
            u.startContainer("samples", "WIG tracks for samples", "WIGsamples", "bigWig")
            for smp in SC.samples:
                u.addWig(smp['wigfile'], smp['name'], smp['name'], "{} atacseq", alwaysZero="on", smoothingWindow="4", visibility="dense")
            u.endContainer()
            u.startContainer("conditions", "WIG tracks for conditions", "WIGconditions", "bigWig")
            for cond in SC.conditions:
                u.addWig(cond['wigfile'], cond['name'], cond['name'], "{} atacseq", alwaysZero="on", smoothingWindow="4", visibility="dense")
            u.endContainer()
            u.startContainer("macspu", "bedGraph tracks for MACS pileup", "BDGpeaks", "bigWig")
            for contr in SC.contrasts:
                u.addBedGraph(contr['macsPileup'], contr['name'] + "-pileup", contr['name'] + "-pileup", "{} - atacseq")
            u.endContainer()
            u.startContainer("macspk", "bedGraph tracks for MACS peaks", "BDGpeaks", "bigWig")
            for contr in SC.contrasts:
                u.addBedGraph(contr['macsNarrow'], contr['name'] + "-peaks", contr['name'] + "-peaks", "{} - atacseq")
            u.endContainer()
            u.startContainer("macssm", "bedGraph tracks for summits", "BDGsummits", "bigWig")
            for contr in SC.contrasts:
                u.addBedGraph(contr['macsSummits'], contr['name'] + "-summits", contr['name'] + "-summits", "{} - atacseq")
            u.endContainer()

        elif self.kind == 'nucleoatac':
            natacdir = "natac/"
            wigs = ["ins" "nucleoatac_signal" "nucleoatac_signal.smooth" "occ"]
            u.startContainer("samples", "WIG tracks for ATAC datasets", "WIGsamples", "bigWig")
            for smp in SC.samples:
                for w in wigs:
                    bwname = natacdir + "/" + smp["name"] + "." + w + ".bw"
                    u.addBigWig(bwname, smp["name"] + "-" + w, smp["name"] + "-" + w, "{} - NucleoATAC")
            u.endContainer()

        # Add hub to Zip file
        ACT.addToZipFile("{}/{}/*".format(ACT.Name, self.dirname))
        # Store hub.txt location
        self.hubpath = self.dirname + "/hub.txt"
        return True

    def Report(self):
        ACT = self.actor
        ACT.scene("UCSC hub")
        
        url = ACT.getConf("url", "Hub")
        if url:
            linkUrl = "http://genome.ucsc.edu/cgi-bin/hgTracks?db={}&hubUrl={}/{}/{}".format(ACT.getConf("genome", "Hub"), url, ACT.Name, self.hubpath)
            ACT.reportf("""Use the following link to display the results in the {}.""", ACT.linkify(linkUrl, "UCSC Genome Browser", target="ucsc"))
        else:
            ACT.reportf("""Use the following link to display the results in the {}. Copy the link location and paste it into the <b>My Hubs</b> form in {}.""",
                        ACT.linkify(self.hubpath, "UCSC Genome Browser", target="ucsc"), 
                        ACT.linkify("http://genome.ucsc.edu/cgi-bin/hgHubConnect", "this page"))
        return True

### For Son's SRA pipeline

class FastqDump(Line):
    name = "fastqdump"

    def Execute(self):
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log

        nsra = 0
        for smp in SC.samples:
            sra = smp['sra']
            outdir = smp['name'] + ".fastqs/"
            smp['fastqdir'] = outdir
            ACT.mkdir(outdir)
            ACT.submit("generic.qsub fastq-dump --split-3 -v --gzip -O {} {} module:sra".format(outdir, sra), done="sra.@.done")
            nsra += 1
        ACT.wait(("sra.@.done", nsra))
        for smp in SC.samples:
            files = glob.glob(smp['fastqdir'] + "*.fastq*")
            # assign left and right fastq files

        return True

REGISTRY = {'test':       testLine,
            'samples':    Samples,
            'rnasamples': RNAseqSamples,
            'trim':       Trimmer,
            'fastqcount': FASTQcounter,
            'bowtie':     Bowtie2,
            'tophat':     TopHat,
            'cufflinks':  Cufflinks,
            'cuffdiff':   Cuffdiff,
            'star':       StarAligner,
            'startx':     StarAlignerTx,
            'bamcount':   BAMcounter,
            'markdup':    Markdup,
            'merge':      BAMmerger,
            'bamcat':     BAMconcatenator,
            'bamcov':     BAMcoverage,
            'kallisto':   Kallisto,
            'rsemquant':  RSEMquant,
            'rsemdiff':   RSEMdiff,
            'bamtowig':   BAMtoWig,
            'macs':       MACScallpeak,
            'tags':       HomerTags,
            'peaks':      HomerPeaks,
            'motifs':     HomerMotifs,
            'diffpeaks':  HomerDiffPeaks,
            'insertsize': InsertSize,
            'csfilter':   CSfilter,
            'mmap':       MMAP,
            'cscall':     CScall,
            'mcomp':      MCOMP,
            'dmr':        DMR,
            'genemeth':   GeneMeth,
            'methplots':  MethPlots,
            'snpgatk':    SNPcallGATK,
            'snpdiff':    SNPdiff,
            'atacqc':     atacQC,
            'hub':        UCSCHub,
}

