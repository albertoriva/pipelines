# ChipSeqActor

# (c) 2016, A. Riva, DiBiG
# ChIPseq processing pipeline, v2.0

import os.path

import Logger
from SampleCollection import SampleCollection

## Initialization

ACT.loadConfiguration(ACT.Arguments[0])
ACT.sc = SampleCollection(ACT.Conf)
SC = ACT.sc
SC.showSamples()

## Get name of STAR index for each sample

defaultstaridx = ACT.getConf("staridx")
for s in SC.samples:
    s['staridx'] = ACT.getConf("staridx", s['name'], default=defaultstaridx)

## Fix paths to fastq files in each readset and store base name.
for rs in SC.readsets:
    rs['left'] = SC.checkPath(rs['left'])
    rs['lbase'] = ACT.setFileExt(os.path.split(rs['left'])[1], "", remove=[".fastq", ".fq", ".gz"]) # basename of left-side reads
    if rs['paired']:
        rs['right'] = SC.checkPath(rs['right'])
        rs['rbase'] = ACT.setFileExt(os.path.split(rs['right'])[1], "", remove=[".fastq", ".fq", ".gz"]) # basename of right-side reads

ACT.title = ACT.getConf("title")
ACT.reference = SC.checkPath(ACT.getConf("reference"))
ACT.btidx = ACT.getConf("btidx") # *** Should verify that this is a bowtie2 index
ACT.setSteps(ACT.getConf("steps"))
ACT.Prefix = ACT.getConf("prefix")

## A few additional things we need
ACT.beforeCounts = "beforeCounts.csv" # Number of aligned reads before dedup
ACT.afterCounts = "afterCounts.csv"   # Number of aligned reads after dedup

## Action starts here
ACT.script(ACT.title, "ChIPseq Analysis Pipeline", "ChIPseq")
ACT.begin(timestamp=False)

LOG = Logger.Logger("chipseq.log", echo='stdout')
LOG.logStart("chipseq")

## Ensure we don't have old .done files lying around
ACT.shell("rm -f *.done *.IN.*")

## Initialize .files
ACT.shell('rm -f .files; echo "index.html *.png *.pdf *.xlsx *.cov.csv" > .files')

## Call trimmomatic on each pair, and then fastqc on the trimmed files
LOG.log("Trimming fastqc files with trimmomatic")
ACT.runTrimmomatic(run=ACT.step("trim"), fastqc=ACT.step("fastqc"))

## Align all reads, and index resulting BAM files
if ACT.step("bowtie"):
    LOG.log("Aligning fastqc files with bowtie")
    ACT.runBowtie(run=ACT.step("bowtie"))
elif ACT.step("star"):
    LOG.log("Aligning fastqc files with STAR")
    ACT.runSTAR(run=ACT.step("star"))

## Remove duplicate reads
LOG.log("Removing duplicate reads with picard")
ACT.doGATK(run=(ACT.step("bowtie") or ACT.step("star")))

## Count number of aligned reads and generate coverage stats
if ACT.step("coverage"):
    ncov = ACT.getBAMcoverage()
    LOG.log("Started {} coverage jobs", ncov)
else:
    ncov = 0

## Merge BAM files for all replicates of each sample
ACT.mergeBAMs()

## Convert BAM files to wig
nwig = ACT.bamToWig()
LOG.log("Started {} bamtowig jobs", nwig)

## Create tag directories for Homer
nhomer = 0
for c in SC.conditions:
    cname = c['name']
    c['tags'] = cname + ".tags.d/"
    c['itags'] = cname + ".itags.d/"
    bams = SC.conditionBAMs(cname)
    LOG.log("BAMs for condition {}: {}", cname, bams)
    ACT.submit("homer.qsub makeTagDirectory {} {}".format(c['tags'], " ".join(bams)), done="homer.@.done")
    nhomer += 1
    ibams = SC.conditionBAMs(cname, role='input')
    LOG.log("BAMs for condition {} (input): {}", cname, ibams)
    ACT.submit("homer.qsub makeTagDirectory {} {}".format(c['itags'], " ".join(ibams)), done="homer.@.done")
    nhomer += 1
ACT.wait(("homer.@.done", nhomer))

## Run peak finding, etc
nhomer = 0
for c in SC.conditions:
    ACT.submit("homer.qsub findPeaks factor {} {}".format(c['tags'], c['itags']), done="homer.@.done")
    nhomer += 1
    ACT.submit("homer.qsub findPeaks histone {} {}".format(c['tags'], c['itags']), done="homer.@.done")
    nhomer += 1
    ACT.submit("homer.qsub findPeaks super {} {}".format(c['tags'], c['itags']), done="homer.@.done")
    nhomer += 1
ACT.wait(("homer.@.done", nhomer))

ACT.wait(("cov.@.done", ncov), ("btw.@.done", nwig))

print "Terminated."
raw_input()


