# An actor that includes a SampleCollection, to process
# multiple samples possibly having multiple replicates each.

import sys
import os.path
from Logger import Logger

# Utils

def cufflinks(outdir, bam, gtf, mask=False):
    cmd = "cuff.qsub cufflinks {} {} {}".format(outdir, bam, gtf)
    if mask:
        cmd += " mask=" + mask
    return cmd

def cuffmerge(outdir, listfile, reference=False, gtf=False):
    cmd = "cuff.qsub cuffmerge {} {}".format(outdir, listfile)
    if reference:
        cmd += " ref=" + reference
    if gtf:
        cmd += " gtf=" + gtf
    return cmd

def cuffdiff(outdir, gtf, labels, bam1, bam2, mask=False):
    cmd = "cuff.qsub cuffdiff {} {} {} {} {}".format(outdir, gtf, ",".join(labels), bam1, bam2)
    if mask:
        cmd += " mask=" + mask
    return cmd

# Main class

from DibigActor import DibigActor

class MultiSampleActor(DibigActor):
    sc = None                   # SampleCollection
    log = Logger(None)          # To avoid errors in scripts that don't explicitly create one

    def runSickle(self, run=True, fastqc=True):
        """Run sickle on the fastq files in all readsets. Pre-trim files are saved under the 'pretrimleft' and 'pretrimright' (for paired-end)
or 'pretrimfastq' (single-end) properties. If `fastqc' is True, run fastqc on the trimmed fastqs."""
        nsickle = 0
        nfastqc = 0
        if fastqc:
            self.mkdir("fastqc")
        for rs in self.sc:
            rs['bamfile'] = self.exclude(rs['name'] + ".bam")
            if rs['paired']:
                fastq1 = rs['left']
                fastq2 = rs['right']

                in1 = os.path.basename(fastq1)
                in2 = os.path.basename(fastq2)
                rs['pretrimleft']  = fastq1
                rs['pretrimright'] = fastq2
                rs['left']  = self.setFileExt(in1, ".sickle.fastq.gz", remove=[".fastq", ".gz"])
                rs['right'] = self.setFileExt(in2, ".sickle.fastq.gz", remove=[".fastq", ".gz"])
                if run:
                    this = self.submit("sickle.qsub pe {} {} {} {} len=50".format(fastq1, fastq2, rs['left'], rs['right']), done="sickle.@.done")
                    nsickle += 1

                rs["fqc1"]   = self.setFileExt(in1, "_fastqc.html", remove=[".fastq", ".gz"])
                rs["fqc2"]   = self.setFileExt(in2, "_fastqc.html", remove=[".fastq", ".gz"])
                if run and fastqc:
                    self.submit("fastqc.qsub {} {}".format(rs["left"], "fastqc"), after=this, done="fqc2.@.done")
                    self.submit("fastqc.qsub {} {}".format(rs["right"], "fastqc"), after=this, done="fqc2.@.done")
                    nfastqc += 2
            else:
                fastq1 = rs['left']
                in1 = os.path.basename(fastq1)
                rs['pretrimleft']  = rs['left']
                rs['left']  = self.setFileExt(in1, ".sickle.fastq.gz", remove=[".fastq", ".gz"])
                if run:
                    this = self.submit("sickle.qsub se {} {} len=50".format(fastq1, rs['left']), done="sickle.@.done")
                    nsickle += 1

                rs["fqc1"]   = self.setFileExt(in1, "_fastqc.html", remove=[".fastq", ".gz"])
                if run and fastqc:
                    self.submit("fastqc.qsub {} {}".format(rs["sickle1"], "fastqc"), after=this, done="fqc2.@.done")
                    nfastqc += 1
        self.wait([("sickle.@.done", nsickle), ("fqc2.@.done", nfastqc)])
            
    def runTrimmomatic(self, run=True, fastqc=True):
        """Run trimmomatic on the fastq files in all readsets. Pre-trim files are saved under the 'pretrimleft' and 'pretrimright' (for paired-end)
or 'pretrimfastq' (single-end) properties. If `fastqc' is True, run fastqc on the trimmed fastqs."""
        ntrimm = 0
        nfastqc = 0
        if fastqc:
            self.mkdir("fastqc")
        for rs in self.sc:
            rs['bamfile'] = self.exclude(rs['name'] + ".bam")
            if rs['paired']:
                fastq1 = rs['left']
                fastq2 = rs['right']

                in1 = os.path.basename(fastq1)
                in2 = os.path.basename(fastq2)
                rs['pretrimleft']  = fastq1
                rs['pretrimright'] = fastq2
                rs['left']  = self.setFileExt(in1, ".trimm.fastq.gz", remove=[".fastq", ".gz"])
                rs['right'] = self.setFileExt(in2, ".trimm.fastq.gz", remove=[".fastq", ".gz"])
                if run:
                    this = self.submit("trimmomatic.qsub {} {} {} {} {} {}".format(
                            fastq1, fastq2, rs['left'], rs['name'] + ".lun.fastq.gz", rs['right'], rs['name'] + ".run.fastq.gz"), 
                                       done="trimm.@.done")
                    ntrimm += 1

                rs["fqc1"]   = self.setFileExt(in1, ".trimm_fastqc.html", remove=[".fastq", ".gz"])
                rs["fqc2"]   = self.setFileExt(in2, ".trimm_fastqc.html", remove=[".fastq", ".gz"])
                if run and fastqc:
                    self.submit("fastqc.qsub {} {}".format(rs["left"], "fastqc"), after=this, done="fqc.@.done")
                    self.submit("fastqc.qsub {} {}".format(rs["right"], "fastqc"), after=this, done="fqc.@.done")
                    nfastqc += 2
            else:
                fastq1 = rs['left']
                in1 = os.path.basename(fastq1)
                rs['pretrimleft']  = rs['left']
                rs['left']  = self.setFileExt(in1, ".trimm.fastq.gz", remove=[".fastq", ".gz"])
                if run:
                    this = self.submit("trimmomatic-se.qsub {} {}".format(fastq1, rs['left']), done="trimm.@.done")
                    ntrimm += 1

                rs["fqc1"]   = self.setFileExt(in1, "_fastqc.html", remove=[".fastq", ".gz"])
                if run and fastqc:
                    self.submit("fastqc.qsub {} {}".format(rs["left"], "fastqc"), after=this, done="fqc.@.done")
                    nfastqc += 1
        self.wait([("trimm.@.done", ntrimm), ("fqc.@.done", nfastqc)])

    def runFastQC(self, run=True, outdir="fastqc"):
        nfastqc = 0
        self.mkdir(outdir)
        for rs in self.sc:
            if rs['paired']:
                fastq1 = rs['left']
                fastq2 = rs['right']

                in1 = os.path.basename(fastq1)
                in2 = os.path.basename(fastq2)

                rs["fqc1"]   = self.setFileExt(in1, "_fastqc.html", remove=[".fastq", ".gz"])
                rs["fqc2"]   = self.setFileExt(in2, "_fastqc.html", remove=[".fastq", ".gz"])
                print "{} -> {}".format(fastq1, rs["fqc1"])
                print "{} -> {}".format(fastq2, rs["fqc2"])
                if run:
                    self.submit("fastqc.qsub {} {}".format(rs["left"], outdir), done="fqc.@.done")
                    self.submit("fastqc.qsub {} {}".format(rs["right"], outdir), done="fqc.@.done")
                    nfastqc += 2
            else:
                fastq1 = rs['left']
                in1 = os.path.basename(fastq1)

                rs["fqc1"]   = self.setFileExt(in1, "_fastqc.html", remove=[".fastq", ".gz"])
                if run:
                    self.submit("fastqc.qsub {} {}".format(rs["left"], outdir), done="fqc.@.done")
                    nfastqc += 1
        self.wait(("fqc.@.done", nfastqc))

    def runBowtie(self, run=True):
        """Run bowtie2 on all the readsets, in single- or paired-end mode. The index is specified in the btidx slot."""
        nbowtie = 0
        for rs in self.sc:
            rs['bamfile'] = self.exclude(rs['name'] + ".bam")
            if run:
                if rs['paired']:
                    self.submit("bowtie2-pe.qsub {} {} {} {}".format(self.btidx, rs['left'], rs['right'], rs['bamfile']), done="bowtie2.@.done")
                else:
                    self.submit("bowtie2.qsub {} {} {}".format(self.btidx, rs['left'], rs['bamfile']), done="bowtie2.@.done")
                nbowtie += 1
        self.wait(("bowtie2.@.done", nbowtie))

    def runSTAR(self, run=True):
        nstar = 0
        for smp in self.sc.samples:
            # Create an output directory for each readset
            for rs in smp['readsets']:
                outdir = rs["name"]
                self.mkdir(outdir)
                rs["bamfile"] = rs["name"] + "/Aligned.sortedByCoord.out.bam" # Name of BAM file produced by STAR
                if run:
                    if rs['paired']:
                        self.submit("star.qsub {} {} {} {}".format(smp["staridx"], rs["left"], rs["right"], outdir), done="star.@.done")
                    else:
                        self.submit("star.qsub {} {} {}".format(smp["staridx"], rs["left"], outdir), done="star.@.done")
                    nstar += 1
        self.wait(("star.@.done", nstar))

    def mergeBAMs(self, run=True, condition=True):
        """If a sample contains multiple replicates, merge its BAM files into a single one. Sets the 'bamfile' entry in the sample to the merged BAM file, or the existing one in case of no replicates. If `condition' is True, also merges the BAM files for all samples in each condition into a single BAM file."""
        nmerge = 0
        for smp in self.sc.samples:
            if len(smp['readsets']) == 1:
                rs = smp['readsets'][0]
                smp['bamfile'] = rs['bamfile']
            else:
                smp['bamfile'] = smp['name'] + ".bam"
                # Only do the merge if the merged file is missing or old
                if self.missingOrStale(smp['bamfile'], other=smp['readsets'][0]['bamfile']):
                    cmdline = "generic.qsub bamtools merge -out " + smp['bamfile']
                    for rs in smp['readsets']:
                        cmdline += " -in " + rs['bamfile']
                    cmdline += " module:bamtools"
                    if run:
                        self.submit(cmdline, done="merge.@.done")
                        nmerge += 1
        self.wait(("merge.@.done", nmerge))
        
        ## Now merge BAMs for all samples in each condition
        nmerge = 0
        for cond in self.sc.conditions:
            condbams = self.getBAMorBAMs(cond['name']).split(",")
            if len(condbams) == 1:
                cond['bamfile'] = condbams[0]
            else:
                cond['bamfile'] = cond['name'] + ".bam"
                if self.missingOrStale(cond['bamfile'], other=condbams[0]):
                    cmdline = "generic.qsub bamtools merge -out " + cond['bamfile']
                    for sb in condbams:
                        cmdline += " -in " + sb
                    cmdline += " module:bamtools"
                    if run:
                        self.submit(cmdline, done="merge.@.done")
                        nmerge += 1
        self.wait(("merge.@.done", nmerge))

    def doGATK(self, run=True):
        ngatk = 0
        for rs in self.sc:
            name = rs['name']
            bamfile = rs['bamfile']
            rs['picard'] = self.exclude(self.setFileExt(bamfile, ".picard.bam"))
            rs['dupmark'] = self.exclude(self.setFileExt(bamfile, ".dupmark.bam"))
            rs['metrics'] = name + ".metrics.txt"
            rs['realigned'] = self.exclude(name + ".realigned.bam")
            rs['bamfile'] = name + ".gatk.bam"

            if run:
                ## bamfile -> addRG
                job1 = self.submit("picard.qsub AddOrReplaceReadGroups I={} O={} RGID=ID_{} RGLB=LB_{} RGPL=ILLUMINA RGPU=PU_{} RGSM=SM_{}".format(
                        bamfile, rs['picard'], name, name, name, name))
                ## addRG -> dupmark
                job2 = self.submit("picard.qsub MarkDuplicates I={} O={} M={} REMOVE_DUPLICATES=true ASSUME_SORTED=true".format(
                        rs['picard'], rs['bamfile'], rs['metrics']), after=job1, done="picard.@.done")
                ## dupmark -> realigned
                # job3 = self.submit("gatk-realigner.qsub {} {} {}".format(self.reference, rs['dupmark'], rs['realigned']), after=job2)
                ## realigned -> bamfile
                #self.submit("picard.qsub FixMateInformation I={} O={} SORT_ORDER=coordinate ASSUME_SORTED=true".format(
                #rs['dupmark'], rs['bamfile']), after=job2, done="picard.@.done")
                ngatk += 1
        self.wait(("picard.@.done", ngatk))
        # Delete unnecessary bam files
        for rs in self.sc:
            self.shell("rm -f {} {} {}".format(rs['picard'], rs['dupmark'], rs['realigned']))

    def postGATK(self, run=True):
        """Index all BAM files."""
        nindex = 0
        for rs in self.sc:
            if run:
                self.submit("bam-index.qsub {}".format(rs['bamfile']), done="index.@.qsub")
                nindex += 1
            # n = self.shell("module load bamtools; bamtools count -in {}".format(rs['bamfile']))
            # if n == '':
            #     sys.stderr.write("Error after GATK step: BAM file `{}' is empty.".format(rs['bamfile']))
            #     rs['bad'] = True
            # else:
            #     rs['fixedReads'] = n
        self.wait(("index.@.qsub", nindex))

    def getBAMcoverage(ACT):
        """Call bamtools coverage on each bam file, followed by chromCoverage.py.
NOTE: this function does not wait for the submitted scripts to run, since they can
take a long time. The calling script should performa a ACT.wait(('cov.@.done', N)) at
the end, where N is the value returned by this function."""
        ncov = 0
        for rs in ACT.sc.readsets:
            bam = "{}.bam".format(rs['name'])
            cov = "{}.cov.csv".format(rs['name'])
            if ACT.missingOrStale(cov, other=bam):
                ncov += 1
                ACT.submit("bamcoverage.qsub {} {} Y".format(bam, cov), done="cov.@.done")
        return ncov

    def countBAMs(self, countsfile, bamfiles, key, slurm=False):
        """If `countsfile' does not exist, write a qsub script to generate it
    by calling bamtools count for each file in `bamfiles'. After it runs, `countsfile'
    will be a tab-delimited file with two columns: readset name, number of reads.
    Then read the contents of `countfile' and set the number of reads as the value
    for the `key' key in each readset."""
        if slurm:
            directives = "#SBATCH --time=10:00:00\n#SBATCH --mem-per-cpu=10G\n#SBATCH --nodes=1\n"
        else:
            directives = "#PBS -l pmem=1gb\n"
        if not os.path.isfile(countsfile):             ## should we re-generate counts?
            with open("countbams.qsub", "w") as out:
                out.write("""#!/bin/bash

{}
    module load bamtools
    rm -f {}
""".format(directives, countsfile))
                for b in bamfiles:
                    out.write("""echo -e "{}\t`bamtools count -in {}`" >> {}\n""".format(b[0], b[1], countsfile))
            self.submit("countbams.qsub", done="countbams.done")
            self.wait("countbams.done")

        ## Now read the counts from the file and assign them to the appropriate readsets
        with open(countsfile, "r") as f:
            for line in f:
                parsed = line.rstrip("\n").split()
                rsname = parsed[0]
                rs = self.sc.findReadset(rsname)
                if rs != None and len(parsed) > 1:
                    rs[key] = int(parsed[1])

    def runCufflinks(self, run=True, mergedir="merge.d/", gtflist="gtflist.txt"):
        """Run cufflinks on all samples. If `mergedir' is specified, also runs cuffmerge on all the GTF files produced by each cufflink."""
        self.cuffmergeOut = mergedir
        self.mkdir(mergedir)

        ncuff = 0
        with open(gtflist, "w") as out:
            for smp in self.sc.samples:
                smp["cuffdir"] = smp["name"] + ".cuff/"
                self.mkdir(smp["cuffdir"])
                out.write("{}/transcripts.gtf\n".format(smp["cuffdir"]))
                if run:
                    self.submit(cufflinks(smp["cuffdir"], smp["bamfile"], smp["cufflinksGTF"], smp["cufflinksMask"]), done="cufflinks.@.done")
                    ncuff += 1
        self.wait(("cufflinks.@.done", ncuff))
        if run:
            self.submit(cuffmerge(mergedir, gtflist, reference=self.reference, gtf=self.cufflinksGTF), done="cuffmerge.done")
            self.wait("cuffmerge.done")

    def getBAMorBAMs(self, what):
        """If `what' is the name of a sample, returns its associated BAM file. If it's the 
name of a condition, returns the BAM files for all its samples as a comma-delimited list."""
        smp = self.sc.findSample(what)
        if smp == None:
            cnd = self.sc.findCondition(what)
            if cnd != None:
                return ",".join([self.sc.findSample(s)['bamfile'] for s in cnd['samples']])
        else:
            return smp['bamfile']

    def runCuffdiff(self, run=True):
        ncont = 0
        for cnt in self.sc.contrasts:
            a = cnt['control']
            b = cnt['test']
            cnt['outdir'] = "{}.vs.{}/".format(b, a)

            bams1 = self.getBAMorBAMs(a)
            bams2 = self.getBAMorBAMs(b)

            if run:
                self.submit(cuffdiff(cnt['outdir'], "merge.d/merged.gtf", [a, b], bams1, bams2), done="cuffdiff.@.done")
                ncont += 1
        self.wait(("cuffdiff.@.done", ncont))

    def runMACS(self, run=True):
        nmacs = 0
        for smp in self.sc.samples:
            smp['bedfile'] = smp['name'] + ".bed"
            if run:
                bams = ",".join([ rs['bamfile'] for rs in smp['readsets'] ])
                self.submit("macs.qsub {} {}".format(smp['name'], bams), done="macs.@.done")
                nmacs += 1
        self.wait(("macs.@.done", nmacs))

    def runHomerPrepare(self, run=True):
        pass
        
