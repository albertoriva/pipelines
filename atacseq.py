# DibigActor

# (c) 2016, A. Riva, DiBiG
# ATACseq processing pipeline, v2.0

from SampleCollection import SampleCollection
from Logger import Logger
from Director import Director

D = Director(ACT)

## Initialization
ACT.loadConfiguration(ACT.Arguments[0])

ACT.staridx = ACT.getConf("staridx")
ACT.bt2idx = ACT.getConf("bt2idx")
ACT.reference = ACT.checkPath(ACT.getConf("reference"))
ACT.cufflinksGTF = ACT.checkPath(ACT.getConf("cufflinksGTF"))
ACT.cufflinksMask = ACT.checkPath(ACT.getConf("cufflinksMask"))
ACT.genesdb = ACT.checkPath(ACT.getConf("genesdb"))
ACT.adapter = ACT.getConf("adapter")
ACT.minlen = ACT.getConf("minlen")
ACT.fdr = ACT.getConf("fdr")
ACT.window = ACT.getConf("window", default="100")
ACT.bamcovd = ACT.getConf("bamcovd", default=5)

## A few additional things we need
ACT.fastqCountsPreTrim = "fastq-counts-pretrim.csv" # Number of input reads before trimming
ACT.fastqCounts = "fastq-counts.csv"                # Number of input reads after trimming
ACT.beforeCounts = "beforeCounts.csv"               # Number of aligned reads before dedup
ACT.afterCounts = "afterCounts.csv"                 # Number of aligned reads after dedup

## Now define the pipeline
D.setSteps(ACT.getConf("steps"))
D.step('samples')
D.step('fastqcount.1', outfile=ACT.fastqCountsPreTrim, propname="fastqCountsPre", delay=True)
D.step('trim', adapter=ACT.adapter, minlen=ACT.minlen)
D.step('fastqcount.2', outfile=ACT.fastqCounts, propname="fastqCounts", delay=True)
D.step('bowtie', extra="optX=2000")
D.step('star')
D.step('bamcount.1', outfile=ACT.beforeCounts, keyname='beforeCounts')
D.step('markdup')
D.step('bamcount.2', outfile=ACT.afterCounts, keyname='afterCounts', delay=True)
D.step('merge', byCondition=True, indexBAM=True)
D.step('atacqc')
D.step('bamcov')
D.step('macs')
D.step('insertsize')
D.step('bamtowig', countsfile=ACT.afterCounts, window=ACT.window)
D.step('hub', kind='atacseq')

D.startAt(ACT.getConf("startat"))
D.stopAt(ACT.getConf("stopAt"))

D.showSteps()

## Action starts here
ACT.script(ACT.title, "ATACseq - Alignment and peak finding")
ACT.begin(timestamp=False)

ACT.initFiles()

D.RunScript()

ACT.cleanup()
