# MultiSampleActor

# (c) 2015, A. Riva, DiBiG
# Differential Methylation Analysis pipeline, v2.0

from SampleCollection import SampleCollection
from Logger import Logger
from Director import Director
import dmaputils

D = Director(ACT)

## Initialization
ACT.loadConfiguration(ACT.Arguments[0])

ACT.title = ACT.getConf("title")
ACT.Prefix = ACT.getConf("label")
ACT.reference = ACT.checkPath(ACT.getConf("reference"))
ACT.siteindex = ACT.checkPath(ACT.getConf("siteindex"))
ACT.mindepth = int(ACT.getConf("mindepth", default="10"))
ACT.minsamples = int(ACT.getConf("minsamples", default="1"))
ACT.strand = ACT.getConf("strand")
ACT.parallel = ACT.getConf("parallel", default="1")
ACT.rrbs = ACT.getConfBoolean("rrbs", default=False)
ACT.setSteps(ACT.getConf("steps"))

## Filtering parameters
ACT.sites = ACT.getConf("sites", section="filter")
ACT.amax = int(ACT.getConf("amax", section="filter", default="10000"))
ACT.fmax = float(ACT.getConf("fmax", section="filter", default="0.05"))

## Differential methylation params
ACT.pval = float(ACT.getConf("pval", section="diff", default=0.05))
ACT.diff = float(ACT.getConf("diff", section="diff", default=0.3))
ACT.genesdb = ACT.getConf("genesdb", section="diff")
ACT.genes = ACT.getConf("genes", section="diff")
ACT.regions = ACT.getConf("regions", section="diff", default='pbd')
ACT.size = ACT.getConf("size", section="diff", default='2000')
ACT.mode = ACT.getConf("mode", section="diff", default='avg')
ACT.nsites = ACT.getConf("nsites", section="diff", default='1')

## Graphics (these should be read from conf file)
ACT.scatterplot_size = 800

## Some more globals (should be read from config file)
ACT.fastqCountsPreTrim = "fastq-counts-pretrim.csv" # Number of input reads before trimming
ACT.fastqCounts =        "fastq-counts.csv"         # Number of input reads after trimming
ACT.fastqCountsFilter =  "fastq-counts.filter.csv"  # Number of input reads after filtering
ACT.alignedCounts =      "alignedCounts.csv"        # Number of aligned reads after de-duplication

## Now define the pipeline
ACT.setSteps(ACT.getConf("steps"))

if ACT.stepPresent('samples'):
    D.add('samples')
if ACT.stepPresent('fastqcount.1'):
    D.add('fastqcount.1', outfile=ACT.fastqCountsPreTrim, propname="fastqCountsPre", delay=True, dry=ACT.stepDry('fastqcount.1'))
if ACT.stepPresent('trim'):
    D.add('trim', dry=ACT.stepDry('trim'))
if ACT.stepPresent('fastqcount.2'):
    D.add('fastqcount.2', outfile=ACT.fastqCounts, propname="fastqCounts", delay=True, dry=ACT.stepDry('fastqcount.2'))
if ACT.stepPresent('csfilter'):
    D.add('csfilter', dry=ACT.stepDry('csfilter'))
if ACT.stepPresent('fastqcount.3'):
    D.add('fastqcount.3', outfile=ACT.fastqCountsFilter, propname="fastqCountsFilter", delay=True, dry=ACT.stepDry('fastqcount.3'))
if ACT.stepPresent('mmap'):
    D.add('mmap', dry=ACT.stepDry('mmap'))
if ACT.stepPresent('markdup'):
    D.add('markdup', dry=ACT.stepDry('markdup'))
if ACT.stepPresent('bamcount'):
    D.add('bamcount', outfile=ACT.alignedCounts, propname='beforeCounts', delay=True, dry=ACT.stepDry('bamcount'))
if ACT.stepPresent('cscall'):
    D.add('cscall', dry=ACT.stepDry('cscall'))
if ACT.stepPresent('mcomp'):
    D.add('mcomp', dry=ACT.stepDry('mcomp'))
if ACT.stepPresent('genemeth'):
    D.add('genemeth', dry=ACT.stepDry('genemeth'))
if ACT.stepPresent('methplots'):
    D.add('methplots', dry=ACT.stepDry('methplots'))

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

## Ensure we don't have old .done files lying around
ACT.shell("rm -f *.done *.IN.*")

## Initialize .files
ACT.shell('rm -f .files; echo -e "*.html\n*.js\n*.css\n*.png\n*.pdf\n*.xlsx\n*.csv\n*.bed\n" > .files')

## Logging
ACT.log.setLogfile(ACT.getConf("logfile"))
ACT.log.setEcho('stdout')
ACT.log.logStart("dmap")

D.RunScript()

ACT.log.logEnd()
ACT.log.close()
