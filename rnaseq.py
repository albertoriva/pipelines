# MultiSampleActor

# (c) 2015, A. Riva, DiBiG
# RNAseq processing pipeline, v2.0

from SampleCollection import SampleCollection
from Logger import Logger
from Director import Director

D = Director(ACT)

## Initialization
ACT.loadConfiguration(ACT.Arguments[0])

ACT.title = ACT.getConf("title")
ACT.Prefix = ACT.getConf("label")
ACT.staridx = ACT.getConf("staridx")
ACT.rsemidx = ACT.getConf("rsemidx")
ACT.reference = ACT.checkPath(ACT.getConf("reference"))
ACT.gtf = ACT.checkPath(ACT.getConf("gtf"))
ACT.cufflinksMask = ACT.checkPath(ACT.getConf("cufflinksMask"))
ACT.fdr = ACT.getConf("fdr")
ACT.fc = float(ACT.getConf("fc"))
ACT.erccdb = ACT.getConf("ERCCdb")

## A few additional things we need
ACT.fastqCountsPreTrim = "fastq-counts-pretrim.csv" # Number of input reads before trimming
ACT.fastqCounts = "fastq-counts.csv"                # Number of input reads after trimming
ACT.beforeCounts = "beforeCounts.csv" # Number of aligned reads before dedup
ACT.afterCounts = "afterCounts.csv"   # Number of aligned reads after dedup
ACT.genomeCounts = "genomeCounts.csv" # Number of aligned reads in genome BAM files

## Now define the pipeline
ACT.setSteps(ACT.getConf("steps"))
print ACT.Steps
if ACT.stepPresent('rnasamples'):
    D.add('rnasamples')
if ACT.stepPresent('fastqcount.1'):
    D.add('fastqcount.1', outfile=ACT.fastqCountsPreTrim, propname="fastqCountsPre", delay=True, dry=ACT.stepDry('fastqcount.1'))
if ACT.stepPresent('trim'):
    D.add('trim', dry=ACT.stepDry('trim'))
if ACT.stepPresent('fastqcount.2'):
    D.add('fastqcount.2', outfile=ACT.fastqCounts, propname="fastqCounts", delay=True, dry=ACT.stepDry('fastqcount.2'))
if ACT.stepPresent('startx'):
    D.add('startx', dry=ACT.stepDry('startx'))
#D.add('merge', byCondition=False, dry=False)
if ACT.stepPresent('bamcount.1'):
    D.add('bamcount.1', outfile=ACT.afterCounts, propname='afterCounts', delay=True, dry=ACT.stepDry('bamcount.1'))
if ACT.stepPresent('bamcount.2'):
    D.add('bamcount.2', outfile=ACT.genomeCounts, propname='genomeCounts', source='genomebam', delay=True, dry=ACT.stepDry('bamcount.2'))
#D.add('bamcat', byCondition=True, dry=False)
if ACT.stepPresent('rsemquant'):
    D.add('rsemquant', dry=ACT.stepDry('rsemquant'))
if ACT.stepPresent('rsemdiff'):
    D.add('rsemdiff', dry=ACT.stepDry('rsemdiff'))
if ACT.stepPresent('bamtowig'):
    D.add('bamtowig', countsfile=ACT.genomeCounts, dry=ACT.stepDry('bamtowig'))

D.startAt(ACT.getConf("startat"))
D.stopAt(ACT.getConf("stopAt"))

D.showSteps()

## Action starts here
ACT.script(ACT.title, "RNAseq - Alignment and differential expression analysis", "RNAseq")
ACT.begin(timestamp=False, copyConf=True)

ACT.log.setLogfile(ACT.getConf("logfile"))
ACT.log.setEcho('stdout')
ACT.log.logStart("rnaseq")

## Ensure we don't have old .done files lying around
ACT.shell("rm -f *.done *.IN.*")

## Initialize .files
ACT.shell('rm -f .files; echo "*.html\n*.png\n*.pdf\n*.xlsx\n*.csv\n*.css\n*.js" > .files')

D.RunScript()

ACT.log.logEnd()
ACT.log.close()
