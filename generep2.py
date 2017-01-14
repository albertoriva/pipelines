# Actor

import os
import sys
import os.path
import genereplib
import Logger

### ToDo:
###   allow using organism name instead of TF file for known ones
###     (e.g., human, mouse)
### Allow specification of mem limits in conf file

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

def filenameNoExt(s):
    return os.path.splitext(os.path.basename(s))[0]

### Script code starts here

## Initialization
ACT.loadConfiguration(ACT.Arguments[0])

ACT.title = ACT.getConf("title")
ACT.Prefix = ACT.getConf("label")
ACT.expressions = ACT.getConf("expressions") # checkPath(ACT.getConf("expressions"))
ACT.tfs = ACT.getConf("tfs")
ACT.translation = ACT.getConf("translation")
ACT.nshuffled = ensureInt(ACT.getConf("nshuffled"))
ACT.samplingsize = ensureInt(ACT.getConf("samplingsize", section="Bootstrap"))
ACT.nrounds = ensureInt(ACT.getConf("nrounds", section="Bootstrap"))
ACT.dpi = ensureFloat(ACT.getConf("dpi", section="Aracne", default="0.1"))
ACT.pvalue = ensureFloat(ACT.getConf("pvalue", section="Aracne", default="1e-7"))

ACT.setSteps(ACT.getConf("steps"))

ACT.script(ACT.title, "GeNeReP - Genetic Network Reconstruction Pipeline", "generep")
ACT.begin(timestamp=False)

LOG = Logger.Logger("generep.log", echo='stdout') # this should be read from conf file
ACT.log = LOG
LOG.logStart("generep")

## Create objects and directories for real and shuffled datasets

ACT.datasets = [ {'name': filenameNoExt(ACT.expressions), 'real': True, 'dirname': 'real'} ]
ACT.real = ACT.datasets[0]
ACT.mkdir("real")

## Remove lines containing too many zeros from expressions file

ACT.shell("cp ../{} real/".format(ACT.expressions))
LOG.log("Expressions file is {}", ACT.expressions)

## Generate randomized datasets

for i in range(ACT.nshuffled):
    name = "shuffled{}".format(i+1)
    dirname = "sh{}".format(i+1)
    ACT.mkdir(dirname)
    ACT.datasets.append({'name': name, 'real': False, 'dirname': dirname})
    LOG.log("Generating randomized dataset {}", name)
    if ACT.step("rand"):
        ACT.shell("module load dibig_tools; cd {}; apple.py random -o {}.txt ../real/{}".format(dirname, name, ACT.expressions))
    
## Perform bootstrap

if ACT.step("bootstrap"):
    nboot = 0
    for d in ACT.datasets:
        LOG.log("Performing bootstrap on file {}.txt, samplingsize={}, rounds={}", d['name'], ACT.samplingsize, ACT.nrounds)
        os.chdir(d['dirname'])
        ACT.submit("generic.qsub apple.py bootstrap -z {} {}.txt {} module:dibig_tools".format(ACT.samplingsize, d['name'], ACT.nrounds), done="../bootstrap.@.done")
        nboot += 1
        os.chdir("..")
    ACT.wait(("bootstrap.@.done", nboot))

## Call aracne

naracne = 0
for d in ACT.datasets:
    LOG.log("Starting ARACNE-AP array for {}", d['name'])
    script = genereplib.submit_aracne_ap_array(d['dirname'] + "/" + d['name'], tfs=ACT.tfs, aend=ACT.nrounds-1, pval=ACT.pvalue, done="../aracne.@.done")
    #print "written for " + d['name']
    #raw_input()
    (scriptdir, scriptfile) = os.path.split(script)
    if ACT.step("aracne"):
        subprocess.call("cd {}; ./{}".format(scriptdir, scriptfile), shell=True)
        naracne += ACT.nrounds
ACT.wait(("aracne.@.done", naracne))

## TODO: compress txt files?

## generate consensus with support=1 using run-step3

ncons = 0
for d in ACT.datasets:
    script = genereplib.write_consensus_script(d['dirname'], d['name'])
    LOG.log("Generating initial consensus with support=1 for {}", d['name'])
    if ACT.step("consensus"):
        os.chdir(d['dirname'])
        ACT.submit(script, done="../consensus.@.done")
        os.chdir("..")
        ncons += 1
ACT.wait(("consensus.@.done", ncons))

## Use compare-mi-histograms on NNN.counts.csv to find optimal support

LOG.log("Comparing real and randomized datasets to determine optimal support")

sum = 0
sumfpr = 0
n = 0
realcounts = ACT.real['dirname'] + "/" + ACT.real['name'] + ".counts.csv"
for d in ACT.datasets:
    if not d['real']:
        outfile = "{}.vs.{}.support.csv".format(ACT.real['name'], d['name'])
        shcounts = d['dirname'] + "/" + d['name'] + ".counts.csv"
        supp = genereplib.compare_mi_histograms(outfile, realcounts, shcounts, maxv=ACT.nrounds+1)
        LOG.log("Comparing {} and {} into {}: maxDiff={} for support={}, FPR={}", realcounts, shcounts, outfile, supp[0], supp[1], supp[2])
        sum += supp[1]
        sumfpr += supp[2]
        n += 1

LOG.log("sumfpr = {}", sumfpr)

optsupport = int(round(1.0 * sum/n))
fprsupport = int(round(1.0 * sumfpr/n))

LOG.log("Re-running consensus generation with optimal support {} (fpr={})", optsupport, fprsupport)

## rerun consensus generation with support found in previous step

ncons = 0
for d in ACT.datasets:
    script = genereplib.write_consensus_script(d['dirname'], d['name'], support=optsupport)
    LOG.log("Generating consensus with support={} for {} into {}.adj", optsupport, d['name'], d['name'])
    if ACT.step("consensus"):
        os.chdir(d['dirname'])
        os.rename(d['name'] + ".adj", d['name'] + ".orig.adj")
        ACT.submit(script, done="../consensus.@.done")
        os.chdir("..")
        ncons += 1
ACT.wait(("consensus.@.done", ncons))

## Use compare_mi_histograms again to find optimal MI threshold

LOG.log("Comparing real and randomized datasets to determine optimal MI")
genereplib.generate_mi_histograms(ACT)

sum = 0
sumfpr = 0
n = 0
realmihist = ACT.real['mihist']
for d in ACT.datasets:
    if not d['real']:
        outfile = "{}.vs.{}.mi.hist.csv".format(ACT.real['name'], d['name'])
        shmihist = d['mihist']
        mi = genereplib.compare_mi_histograms(outfile, realmihist, shmihist)
        LOG.log("Comparing {} and {}: maxDiff={} for mi={}, FPR={}", realmihist, shmihist, mi[0], mi[1], mi[2])
        sum += mi[1]
        sumfpr += mi[2]
        n += 1

optmi = sum/n
fprmi = sumfpr/n

LOG.log("MI threshold: {} (fpr={})", optmi, fprmi)

## Now filter the adj files with this threshold
nfilt = 0
for d in ACT.datasets:
    LOG.log("Filtering dataset {} with MI={}", d['name'], optmi)
    script = genereplib.write_filter_script(d['dirname'], d['name'], mi=optmi)
    os.chdir(d['dirname'])
    ACT.submit(script, done="../filter.@.done")
    os.chdir("..")
    nfilt += 1
ACT.wait(("filter.@.done", nfilt))

## Use compare_mi_histograms again to find optimal sum(MI) threshold

LOG.log("Comparing real and randomized datasets to determine optimal sum(MI)")
genereplib.generate_mi_histograms(ACT, summi=True)

sum = 0
sumfpr = 0
n = 0
realmihist = ACT.real['summihist']
for d in ACT.datasets:
    if not d['real']:
        outfile = "{}.vs.{}.summi.hist.csv".format(ACT.real['name'], d['name'])
        shmihist = d['summihist']
        mi = genereplib.compare_mi_histograms(outfile, realmihist, shmihist)
        LOG.log("Comparing {} and {}: maxDiff={} for sumMI={}, FPR={}", realmihist, shmihist, mi[0], mi[1], mi[2])
        sum += mi[1]
        sumfpr += mi[2]
        n += 1

optsummi = sum/n
fprsummi = sumfpr/n

LOG.log("SumMI threshold: {} (fpr={})", optsummi, fprsummi)

## Now filter the adj files again with this threshold
nfilt = 0
for d in ACT.datasets:
    LOG.log("Filtering dataset {} with sum(MI)={}", d['name'], optsummi)
    script = genereplib.write_filter_script(d['dirname'], d['name'], summi=optsummi)
    os.chdir(d['dirname'])
    ACT.submit(script, done="../filter.@.done")
    os.chdir("..")
    nfilt += 1
ACT.wait(("filter.@.done", nfilt))

## Write stats on all adj files after each filtering step
genereplib.all_adj_stats(ACT, "adj-stats.csv")
with open("adj-stats.csv", "r") as f:
    stats = f.read()
LOG.log("ADJ stats after each filtering step:\n{}", stats)

## Convert final adj file to cytoscape
final = ACT.real['dirname'] + "/" + ACT.real['name'] + ".summi.adj"
cyto = ACT.real['name'] + ".cyto.csv"
LOG.log("Converting final file {} to cytoscape file {}", final, cyto)
ACT.shell("module load dibig_tools; apple.py convert ac {} {}".format(final, cyto))

# If translation file specified, do translation
if ACT.translation:
    cytopre = cyto
    cyto = ACT.real['name'] + ".tr.csv"
    LOG.log("Translating gene names in {} to {} using map {}", cytopre, cyto, ACT.translation)
    ACT.shell("module load dibig_tools; apple.py translate {} {} {}".format(ACT.translation, cytopre, cyto))
    finalpre = final
    final = ACT.real['name'] + ".tr.adj"
    LOG.log("Converting translated file {} to {}", cyto, final)
    ACT.shell("module load dibig_tools; apple.py convert ca {} {}".format(cyto, final))

# Produce connections file and print top genes.
connections = ACT.real['name'] + ".conn.csv"
ACT.shell("apple.py convert co {} {}".format(cyto, connections))
tophubs = ACT.shell("cut -f 1,2 {} | head -50".format(connections))
LOG.log("Top 50 hubs by number of connections:\n" + tophubs)

LOG.logEnd()
LOG.close()

### Report generation

scene = 1
ACT.scene(scene, "General configuration")
ACT.reportf("""Expressions file: <b>{}</b><br>
Number of randomizations: <b>{}</b><br>
Bootstrap rounds: <b>{}</b><br>
Bootstrap sampling size: <b>{}</b><br>
Aracne DPI: <b>{}</b><br>
Aracne P-value: <b>{}</b><br>""".format(ACT.expressions, ACT.nshuffled, ACT.nrounds, ACT.samplingsize, ACT.dpi, ACT.pvalue))

scene += 1
ACT.scene(scene, "Randomization and bootstrap")
ACT.reportf("""<P><b>{}</b> randomized datasets were generated by random shuffling.</P>""".format(ACT.nshuffled))
ACT.reportf("""<P>Each randomized dataset was bootstrapped into <b>{}</b> new files, with a sample size of <b>{}</b>.</p>""".format(ACT.nrounds, ACT.samplingsize))

scene += 1
ACT.scene(scene, "Aracne analysis and consensus generation")
ACT.reportf("""<P>ADJ files were generated using <b>Aracne</b> with the following parameters: DPI={}, P-value={}.</P>""".format(ACT.dpi, ACT.pvalue))
ACT.reportf("""<P>The ADJ files were merged into a consensus network for each sample. The following table lists all consensus networks.</P>""")
rows = []
for d in ACT.datasets:
    rows.append( [ d['name'], ACT.linkify(d['dirname'] + "/" + d['name'] + ".adj", d['name'] + ".adj") ] )
ACT.table(rows, header=['Dataset', 'Consensus'], align="HL")

scene += 1
ACT.scene(scene, "Filtering")
ACT.reportf("""<P>The real and randomized consensus networks were compared to determine the optimal support, MI and sum(MI) threshold. 
The following values were found:</P>""")
ACT.table( [ ["Support", optsupport, "{:.3f}".format(fprsupport)],
             ["MI", "{:.3f}".format(float(optmi)), "{:.3f}".format(fprmi)],
             ["sum(MI)", "{:.3f}".format(float(optsummi)), "{:.3f}".format(fprsummi)] ],
           header = ["Parameter", "Optimal", "FPR"],
           align="HRR")

finalstats = subprocess.check_output("module load dibig_tools; apple.py stats {}".format(final), shell=True)
finalstats = finalstats.split()
ACT.reportf("""<P>Details of consensus network:</P>""")
ACT.table( [ [ "Number of rows:", finalstats[5] ],
             [ "Total number of edges:", finalstats[6] ],
             [ "Average number of edges:", "{:.2f}".format(float(finalstats[7])) ] ],
           align="HR")
ACT.file(final, description="Filtered consensus network (ADJ format).")
ACT.file(cyto, description="Filtered consensus network (Cytoscape format).")
