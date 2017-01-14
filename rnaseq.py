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
ACT.bt2idx = ACT.getConf("bt2idx")
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
D.setSteps(ACT.getConf("steps"))
D.step('rnasamples')
D.step('fastqcount.1', outfile=ACT.fastqCountsPreTrim, propname="fastqCountsPre", delay=True)
D.step('trim')
D.step('fastqcount.2', outfile=ACT.fastqCounts, propname="fastqCounts", delay=True)
D.step('startx')
#D.step('tophat')
#D.step('cufflinks')
#D.step('cuffdiff')
D.step('merge', indexBAM=True)
D.step('bamcount.1', outfile=ACT.afterCounts, propname='afterCounts', delay=True)
D.step('bamcount.2', outfile=ACT.genomeCounts, propname='genomeCounts', source='genomebam', delay=True)
D.step('rsemquant')
D.step('rsemdiff')
D.step('bamtowig', countsfile=ACT.genomeCounts)

D.startAt(ACT.getConf("startat"))
D.stopAt(ACT.getConf("stopAt"))

D.showSteps()

## Action starts here
ACT.script(ACT.title, "RNAseq - Alignment and differential expression analysis", "RNAseq")
ACT.begin(timestamp=False, copyConf=True)

ACT.initFiles()

D.RunScript()

ACT.cleanup()
