# DibigActor

# (c) 2015, A. Riva, DiBiG
# Differential Methylation Analysis pipeline, v2.0

### TODO: convert BED to bedGraph automatically with sliding window
### add conf file parameters for dmr

from SampleCollection import SampleCollection
from Logger import Logger
from Director import Director

D = Director(ACT)

## Initialization
ACT.loadConfiguration(ACT.Arguments[0])

ACT.title = ACT.getConf("title")
ACT.Prefix = ACT.getConf("label")
ACT.reference = ACT.checkPath(ACT.getConf("reference"))
ACT.siteindex = ACT.checkPath(ACT.getConf("siteindex"))
ACT.mindepth = ACT.getConfInt("mindepth", default=10)
ACT.minsamples = ACT.getConfInt("minsamples", default=1)
ACT.strand = ACT.getConf("strand")
ACT.parallel = ACT.getConf("parallel", default="1")
ACT.rrbs = ACT.getConf("RRBS", default=None)
ACT.setSteps(ACT.getConf("steps"))

## Trimming
ACT.adapter = ACT.getConf("adapter")
ACT.minlen = ACT.getConf("minlen")

## Filtering parameters
ACT.sites = ACT.getConf("sites", section="filter")
ACT.amax = ACT.getConfInt("amax", section="filter", default=10000)
ACT.fmax = ACT.getConfFloat("fmax", section="filter", default=0.05)

## Differential methylation params
ACT.pval = ACT.getConfFloat("pval", section="diff", default=0.05)
ACT.diff = ACT.getConfFloat("diff", section="diff", default=0.3)
ACT.genesdb = ACT.getConf("genesdb", section="diff")
ACT.genes = ACT.getConf("genes", section="diff")
ACT.regions = ACT.getConfList("regions", section="diff", default='pbd')
ACT.size = ACT.getConfList("size", section="diff", default='2000')
ACT.mode = ACT.getConfList("mode", section="diff", default='avg')
ACT.nsites = ACT.getConfList("nsites", section="diff", default='1')

## DMR parameters
ACT.dmr_winsize = ACT.getConfInt("winsize", "dmr", 100)
ACT.dmr_mincov  = ACT.getConfInt("mincov", "dmr", 4)
ACT.dmr_minsites = ACT.getConfInt("minsites", "dmr", 4)
ACT.dmr_diff     = ACT.getConfFloat("diff", "dmr", 0.2)
ACT.dmr_pval     = ACT.getConfFloat("pval", "dmr", 0.01)
ACT.dmr_gap      = ACT.getConfInt("gap", "dmr", 1)

## Graphics (these should be read from conf file)
ACT.scatterplot_size = 800

## Some more globals (should be read from config file)
ACT.fastqCountsPreTrim = "fastq-counts-pretrim.csv" # Number of input reads before trimming
ACT.fastqCounts =        "fastq-counts.csv"         # Number of input reads after trimming
ACT.fastqCountsFilter =  "fastq-counts.filter.csv"  # Number of input reads after filtering
ACT.alignedCounts =      "alignedCounts.csv"        # Number of aligned reads after de-duplication

## Now define the pipeline
D.setSteps(ACT.getConf("steps"))

D.step('samples')
D.step('fastqcount.1', outfile=ACT.fastqCountsPreTrim, propname="fastqCountsPre", delay=True)
D.step('trim', adapter=ACT.adapter, minlen=ACT.minlen)
D.step('fastqcount.2', outfile=ACT.fastqCounts, propname="fastqCounts", delay=True)
D.step('csfilter')
D.step('fastqcount.3', outfile=ACT.fastqCountsFilter, propname="fastqCountsFilter", delay=True)
D.step('mmap')
#if ACT.stepPresent('markdup') and ACT.rrbs == None:
D.step('markdup')
D.step('bamcount', outfile=ACT.alignedCounts, propname='beforeCounts', delay=True)
D.step('bamcov')
D.step('cscall')
D.step('mcomp')
D.step('dmr')
D.step('genemeth')
D.step('methplots')

# if ACT.stepPresent('bamcount.1'):
#     D.add('bamcount.1', outfile=ACT.beforeCounts, propname='beforeCounts', delay=False, dry=ACT.stepDry('bamcount.1'))
# #D.add('merge', byCondition=False, dry=False)
# if ACT.stepPresent('bamcount.2'):
#     D.add('bamcount.2', outfile=ACT.afterCounts, propname='afterCounts', delay=True, dry=ACT.stepDry('bamcount.2'))
# if ACT.stepPresent('bamcount.3'):
#     D.add('bamcount.3', outfile=ACT.genomeCounts, propname='genomeCounts', source='genomebam', delay=True, dry=ACT.stepDry('bamcount.3'))
# #D.add('bamcat', byCondition=True, dry=False)
# if ACT.stepPresent('bamtowig'):
#     D.add('bamtowig', countsfile=ACT.genomeCounts, dry=ACT.stepDry('bamtowig'))

D.startAt(ACT.getConf("startat"))
D.stopAt(ACT.getConf("stopAt"))

D.showSteps()

## Action starts here
ACT.script(ACT.title, "DMAP - Differential Methylation Analysis Pipeline", "dmap")
ACT.begin(timestamp=False, copyScript=False, copyConf=True)

ACT.initFiles()

D.RunScript()

ACT.cleanup()
