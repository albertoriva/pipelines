# (c) 2016, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import os
import os.path
import subprocess

import Utils
from Lines import Line
import genereputils

##

class Dataset():
    name = ""
    dirname = ""
    real = True
    
    # Initial adj, counts, and network file
    adj = ""
    counts = ""
    network = ""

    # Adj, counts, and network file after support filtering
    adj_supp = ""
    counts_supp = ""
    network_supp = ""

    # Histograms
    mihist = ""                 # MI histogram
    summihist = ""              # sum(MI) histogram
    finalhist = ""              # final histogram

    def __init__(self, name, dirname, real=True):
        self.name = name
        self.dirname = dirname
        self.real = real

    def path(self, what=None):
        if not what:
            what = self.name
        return self.dirname + "/" + what

    """Write the qsub file to generate the ADJ consensus for the specified `support' (defaulting to 1 if None). If `network' is True, also generates the network file.
"""

    def write_consensus_script(self, outfile="consensus.qsub", support=None, network=False):
        outpath = self.dirname + "/" + outfile
        with open(outpath, "w") as out:
            out.write("""#!/bin/bash
#SBATCH --mem-per-cpu=90G
#SBATCH --time=60:00:00

module load dibig_tools

""")
            if support:
                self.adj_supp     = "{}.supp={}.adj".format(self.name, support)
                self.counts_supp  = "{}.counts.supp={}.csv".format(self.name, support)
                self.network_supp = "{}.netwrk.supp={}.csv".format(self.name, support)
                nopt = "-d " + self.network.supp if network else ""
                out.write("""apple.py consensus -s {} -c {} {} {} {}-*.adj*\n""".format(
                        support, self.counts_supp, nopt, self.adj_supp, self.name))
            else:
                self.adj     = self.name + ".adj"
                self.counts  = "{}.counts.csv".format(self.name)
                self.network = "{}.netwrk.csv".format(self.name)
                nopt = "-d " + self.network if network else ""
                out.write("""apple.py consensus -s 1 -c {} {} {} {}-*.adj*\n""".format(self.counts, nopt, self.adj, self.name))
        return outfile

    def write_filter_script(self, outfile="filter.qsub", mi=None, summi=None, log=None):
        directory = self.dirname
        prefix = self.name
        outpath = self.path(outfile)
        if mi:
            inadj = self.adj_supp                                          # "{}.adj".format(prefix)
            outadj = self.adj_mi = "{}.mi.adj".format(prefix)
            flag = ""
        elif summi:
            inadj = self.adj_mi                                            # "{}.mi.adj".format(prefix)
            outadj = self.adj_summi = "{}.summi.adj".format(prefix)
            flag = "-t"
            mi = summi
        else:
            print "Error: one of `mi' and `summi' must be specified!"
            sys.exit(1)
        if log:
            log.log("Writing {}: {} => {}", outpath, inadj, outadj)
        with open(outpath, "w") as out:
            out.write("""#!/bin/bash
#SBATCH --mem=40G
#SBATCH --time=60:00:00

module load dibig_tools

apple.py filter {} -o {} {} {}
""".format(flag, outadj, inadj, mi))
        return outfile

## Generep

class GenerepInit(Line):
    name = "init"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        ACT.datasets = [ Dataset(Utils.filenameNoExt(ACT.expressions), 'real', True) ]
        ACT.real = ACT.datasets[0]
        ACT.mkdir("real")

        ## (Remove lines containing too many zeros from expressions file)

        #ACT.shell("cp ../{} real/".format(ACT.expressions))
        source = "../" + ACT.expressions
        dest   = "real/" + ACT.expressions
        LOG.log("Expressions file is {}", ACT.expressions)
        if ACT.missingOrStale(dest, other=source):
            (nin, nout) = genereputils.removeZeros(source, dest, genecols=ACT.genecols)
            LOG.log("Removed {} genes with too many missing values, {} genes retained.", nin-nout, nout)
        return True

class GenerepRand(Line):
    name = "rand"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        for i in range(ACT.nshuffled):
            name = "shuffled{}".format(i+1)
            dirname = "sh{}".format(i+1)
            ACT.mkdir(dirname)
            ACT.datasets.append(Dataset(name, dirname, real=False))
            LOG.log("Generating randomized dataset {}", name)
            if not self.dry:
                ACT.shell("module load dibig_tools; cd {}; apple.py random -o {}.txt -gc {} ../real/{}".format(dirname, name, ACT.genecols, ACT.expressions))
        return True

class GenerepBootstrap(Line):
    name = "bootstrap"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        nboot = 0
        for d in ACT.datasets:
            LOG.log("Performing bootstrap on file {}.txt, samplingsize={}, rounds={}", d.name, ACT.samplingsize, ACT.nrounds)
            os.chdir(d.dirname)
            if not self.dry:
                ACT.submit("generic.qsub apple.py bootstrap -z {} -gc {} {}.txt {} module:dibig_tools".format(ACT.samplingsize, ACT.genecols, d.name, ACT.nrounds), done="../bootstrap.@.done")
                nboot += 1
            os.chdir("..")
        return ACT.wait(("bootstrap.@.done", nboot))

class GenerepAracne(Line):
    name = "aracne"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        naracne = 0
        for d in ACT.datasets:
            LOG.log("Starting ARACNE array for {}", d.name)
            script = genereputils.submit_aracne_array(d.path(), tfs=ACT.tfs, aend=ACT.nrounds-1, pval=ACT.pvalue, done="../aracne.@.done", qsuboptions=ACT.qsuboptions)
            (scriptdir, scriptfile) = os.path.split(script)
            if not self.dry:
                subprocess.call("cd {}; ./{}".format(scriptdir, scriptfile), shell=True)
                naracne += ACT.nrounds
        return ACT.wait(("aracne.@.done", naracne))

class GenerepAracneAP(Line):
    name = "aracneap"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        naracne = 0
        for d in ACT.datasets:
            LOG.log("Starting ARACNE-AP array for {}", d['name'])
            script = genereputils.submit_aracne_ap_array(d['dirname'] + "/" + d['name'], tfs=ACT.tfs, aend=ACT.nrounds-1, pval=ACT.pvalue, done="../aracne.@.done")
            (scriptdir, scriptfile) = os.path.split(script)
            if not self.dry:
                subprocess.call("cd {}; ./{}".format(scriptdir, scriptfile), shell=True)
                naracne += ACT.nrounds
        return ACT.wait(("aracne.@.done", naracne))

class GenerepOrigAracne(Line):
    name = "origaracne"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        ACT.norigaracne = 0
        for d in ACT.datasets:
            LOG.log("Starting original ARACNE-AP array for {}", d.name)
            script = genereputils.submit_aracne_ap_array(d.path(), tfs=ACT.tfs, aend=ACT.nrounds-1, pval=ACT.pvalue, done="../origaracne.@.done", orig=True)
            (scriptdir, scriptfile) = os.path.split(script)
            if not self.dry:
                subprocess.call("cd {}; ./{}".format(scriptdir, scriptfile), shell=True)
                ACT.norigaracne += 1
        return True

    def PostExecute(self):
        ACT = self.actor
        return ACT.wait(("origaracne.@.done", ACT.norigaracne))

class GenerepConsensus(Line):
    name = "consensus"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        ncons = 0
        for d in ACT.datasets:
            script = d.write_consensus_script()
            LOG.log("Generating initial consensus with support=1 for {}", d.name)
            LOG.log("Current .adj file: {}", d.adj)
            if not self.dry:
                os.chdir(d.dirname)
                ACT.submit(script, done="../consensus.@.done")
                os.chdir("..")
                ncons += 1
        return ACT.wait(("consensus.@.done", ncons))

### Histogram-based filtering
desc = """
At the beginning of this step we have one consensus adj file in each subdir, called NAME.adj,
and one counts file, called NAME.counts.csv.

ADJ stored in:  d.adj

1. compare_mi_histograms() on NAME.counts.csv for real and each randomized dataset. Results in:
   REAL.vs.RND.support.csv. Determines optimal support (optsupport, fprsupport).

2. Generate consensus using new support S. New adj called NAME.supp=S.adj. Stored in d.adj_supp

3. generate_mi_histograms() with summi=False and final=False. Writes 

"""

class GenerepHistograms(Line):
    name = "histograms"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        filter1 = genereputils.GRFilterSupport(ACT)
        filter2 = genereputils.GRFilterMI(ACT)
        filter3 = genereputils.GRFilterSumMI(ACT)
        filter4 = genereputils.GRFilterFinal(ACT)

        LOG.log("*** Filter 1")
        filter1.run(dry=self.dry)
        LOG.log("*** Filter 2")
        filter2.run(dry=self.dry)
        LOG.log("*** Filter 3")
        filter3.run(dry=self.dry)
        LOG.log("*** Filter 4")
        filter4.run(dry=self.dry)
        return True

    def ExecuteOld(self):
        ACT = self.actor
        LOG = ACT.log

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
                supp = genereputils.compare_mi_histograms(outfile, realcounts, shcounts, maxv=ACT.nrounds+1)
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
            (script, adj) = genereputils.write_consensus_script(d['dirname'], d['name'], support=optsupport)
            LOG.log("Generating consensus with support={} for {}", optsupport, d['name'])
            d['adj'] = adj
            LOG.log("Current .adj file: {}.adj", adj)
            if not self.dry:
                os.chdir(d['dirname'])
                # os.rename(d['name'] + ".adj", d['name'] + ".orig.adj")
                ACT.submit(script, done="../consensus.@.done")
                os.chdir("..")
                ncons += 1
        ACT.wait(("consensus.@.done", ncons))

        ## Use compare_mi_histograms again to find optimal MI threshold

        LOG.log("Comparing real and randomized datasets to determine optimal MI")
        genereputils.generate_mi_histograms(ACT)

        sum = 0
        sumfpr = 0
        n = 0
        realmihist = ACT.real['mihist']
        for d in ACT.datasets:
            if not d['real']:
                outfile = "{}.vs.{}.mi.hist.csv".format(ACT.real['name'], d['name'])
                shmihist = d['mihist']
                mi = genereputils.compare_mi_histograms(outfile, realmihist, shmihist)
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
            script = genereputils.write_filter_script(d['dirname'], d['name'], mi=optmi)
            if not self.dry:
                os.chdir(d['dirname'])
                ACT.submit(script, done="../filter.@.done")
                os.chdir("..")
                nfilt += 1
        ACT.wait(("filter.@.done", nfilt))

        ## Use compare_mi_histograms again to find optimal sum(MI) threshold

        LOG.log("Comparing real and randomized datasets to determine optimal sum(MI)")
        genereputils.generate_mi_histograms(ACT, summi=True)

        sum = 0
        sumfpr = 0
        n = 0
        realmihist = ACT.real['summihist']
        for d in ACT.datasets:
            if not d['real']:
                outfile = "{}.vs.{}.summi.hist.csv".format(ACT.real['name'], d['name'])
                shmihist = d['summihist']
                mi = genereputils.compare_mi_histograms(outfile, realmihist, shmihist)
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
            script = genereputils.write_filter_script(d['dirname'], d['name'], summi=optsummi)
            if not self.dry:
                os.chdir(d['dirname'])
                ACT.submit(script, done="../filter.@.done")
                os.chdir("..")
                nfilt += 1
        ACT.wait(("filter.@.done", nfilt))

        ## Write final MI histogram
        LOG.log("Comparing real and randomized datasets to generate final MI histogram.")
        genereputils.generate_mi_histograms(ACT, final=True)
        ## *** TODO: average final histograms into single one, or get apple.py to average them.
        
## Write stats on all adj files after each filtering step
        genereputils.all_adj_stats(ACT, "adj-stats.csv")
        with open("adj-stats.csv", "r") as f:
            stats = f.read()
        LOG.log("ADJ stats after each filtering step:\n{}", stats)

## Convert final adj file to cytoscape
        final = ACT.real['dirname'] + "/" + ACT.real['name'] + ".summi.adj"
        cyto = ACT.real['name'] + ".cyto.csv"
        LOG.log("Converting final file {} to cytoscape file {}", final, cyto)
        ACT.shell("module load dibig_tools; apple.py convert ac -f {} {}".format(final, cyto))

## If translation file specified, do translation

        if ACT.translation:
            cytopre = cyto
            cyto = ACT.real['name'] + ".tr.csv"
            LOG.log("Translating gene names in {} to {} using map ../{}", cytopre, cyto, ACT.translation)
            ACT.shell("module load dibig_tools; apple.py translate ../{} {} {}".format(ACT.translation, cytopre, cyto))
            finalpre = final
            final = ACT.real['name'] + ".tr.adj"
            LOG.log("Converting translated file {} to {}", cyto, final)
            ACT.shell("module load dibig_tools; apple.py convert ca {} {}".format(cyto, final))

## Produce connections file and CX file, and print top genes.
        connections = ACT.real['name'] + ".conn.csv"
        ACT.shell("module load dibig_tools; apple.py convert co {} {}".format(cyto, connections))
        cxfile = ACT.real['name'] + ".cx"
        
        if len(ACT.attributes) > 0:
            with open("attributes.txt", "w") as out:
                for a in ACT.attributes:
                    out.write("{}: {}\n".format(a[0].capitalize(), a[1]))
            afile = "-a attributes.txt"
        else:
            afile = ""

        ACT.shell("module load dibig_tools; apple.py convert ax {} {} {}".format(afile, final, cxfile))
        tophubs = ACT.shell("cut -f 1,2 {} | head -50".format(connections))
        LOG.log("Top 50 hubs by number of connections:\n" + tophubs)

        return True

# Compression
class GenerepCompress(Line):
    name = "compress"
    nzips = 30                  # Number of files to zip in single gzip.qsub job

    def Execute(self):
        ngz = 0                     # Number of gzip.qsub jobs submitted
        ACT = self.actor
        ACT.shell("rm -f comprfiles.tmp; find . -name \*.adj -size +10M > comprfiles.tmp; find . -name \*.txt -size +10M >> comprfiles.tmp")
        nfiles = int(ACT.fileLines("comprfiles.tmp"))
        for i in range(1, nfiles, self.nzips):
            if not self.dry:
                ACT.submit("gzip.qsub @comprfiles.tmp {} +{}".format(i, self.nzips), done="compr.@.done")
                ngz += 1

        return ACT.wait(("compr.@.done", ngz))

REGISTRY = {'init': GenerepInit,
            'randomize': GenerepRand,
            'bootstrap': GenerepBootstrap,
            'aracne': GenerepAracne,
            'aracneAP': GenerepAracneAP,
            'origaracne': GenerepOrigAracne,
            'consensus': GenerepConsensus,
            'histograms': GenerepHistograms,
            'compress': GenerepCompress}

