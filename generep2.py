# DibigActor

import os
import os.path
from Director import Director

### ToDo:
###   allow using organism name instead of TF file for known ones
###     (e.g., human, mouse)
### Allow specification of mem limits in conf file
### Clean up handling of adj files in histogram step (plus compression)
### Use FPR instead of MI in ax conversion

def ensureInt(s):
    try:
        return int(s)
    except ValueError:
        return None

def ensureFloat(s):
    try:
        return float(s)
    except ValueError:
        return False

D = Director(ACT, "generepLib")

### Script code starts here

## Initialization
ACT.loadConfiguration(ACT.Arguments[0])

ACT.title = ACT.getConf("title")
ACT.Prefix = ACT.getConf("label")
ACT.expressions = ACT.getConf("expressions") # checkPath(ACT.getConf("expressions"))
ACT.tfs = ACT.getConf("tfs")
ACT.translation = ACT.getConf("translation")
ACT.genesdb = ACT.getConf("genesdb")
ACT.genecols = ensureInt(ACT.getConf("genecols", default=2))
ACT.nshuffled = ensureInt(ACT.getConf("nshuffled"))
ACT.samplingsize = ensureInt(ACT.getConf("samplingsize", section="Bootstrap"))
ACT.nrounds = ensureInt(ACT.getConf("nrounds", section="Bootstrap"))
ACT.dpi = ensureFloat(ACT.getConf("dpi", section="Aracne", default="0.1"))
ACT.pvalue = ensureFloat(ACT.getConf("pvalue", section="Aracne", default="1e-7"))
ACT.qsuboptions = ACT.getConf("qsuboptions")
ACT.tophubs = ACT.getConf("tophubs", default=50)

## Options for CX file generation
ACT.attributes = ACT.getConfAll("CX")

## Define the pipeline
D.setSteps(ACT.getConf("steps"))
D.step("init")
D.step("randomize")
D.step("origaracne")
D.step("bootstrap")
D.step("aracne")
#D.step("aracneAP")
D.step("consensus")
D.step("histograms")
D.step("compress")

D.startAt(ACT.getConf("startat"))
D.stopAt(ACT.getConf("stopAt"))

D.showSteps()

## Action starts here
ACT.script(ACT.title, "GeNeReP - Genetic Network Reconstruction Pipeline", "generep")
ACT.begin(timestamp=False)

ACT.initFiles()

D.RunScript()

ACT.cleanup()


