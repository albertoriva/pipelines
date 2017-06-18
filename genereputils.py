import os
import glob
import os.path
import subprocess
from datetime import date

### Functions implementing script steps

def submit_aracne_array(infile, outfile="aracne-all.sh", tfs=False, astart=0, aend=99, amax=40, pval="1e-7", dpi="0.1", done=None, qsuboptions=None):
    (dirname, filename) = os.path.split(infile)
    (prefix, ext) = os.path.splitext(filename)

    outpath = dirname + "/" + outfile

    if done:
        done = "-done " + done
    else:
        done = ""
    if qsuboptions:
        qsuboptions = "-o " + qsuboptions
    else:
        qsuboptions = ""
    print "Writing {}".format(outpath)
    with open(outpath, "w") as out:
        out.write("""#!/bin/bash
module load dibig_tools
""")
        out.write("submit -t {}-{}%{} {} {} aracne-array.qsub {} {} pval={} dpi={}".format(astart, aend, amax, qsuboptions, done, prefix, prefix, pval, dpi))
        if tfs:
            if tfs[0] != '/':
                tfs = "../" + tfs
            out.write(" tfs={}".format(tfs))
        out.write("\n")
    subprocess.call("chmod +x " + outpath, shell=True)
    return outpath

def submit_aracne_ap_array(infile, outfile="aracne-all.sh", tfs=False, astart=0, aend=99, amax=40, pval="1e-7", orig=False, done=None):
    (dirname, filename) = os.path.split(infile)
    (prefix, ext) = os.path.splitext(filename)

    outpath = dirname + "/" + outfile

    if done:
        done = "-done " + done
    else:
        done = ""

    print "Writing {}".format(outpath)
    with open(outpath, "w") as out:
        out.write("""#!/bin/bash
module load dibig_tools
""")
        if tfs:
            if tfs[0] == '/':
                tfarg = " tfs={}".format(tfs)
            else:
                tfarg = " tfs=../{}".format(tfs)
        else:
            tfarg = ""

        if orig:
            prefix = "AP_" + prefix
            for f in glob.glob("{}/*".format(prefix)):
                os.unlink(f)
            out.write("J=`submit -t {}-{}%{} aracne-ap-array.qsub e={} o={} {} mode=boot`\n".format(astart, aend, amax, filename + ".txt", prefix, tfarg))
            out.write("submit {} -after $J aracne-ap-array.qsub o={} adj={} mode=consolidate\n".format(done, prefix, filename + ".txt"))
        else:
            out.write("submit -t {}-{}%{} {} aracne-ap-array.qsub e={} o={} pvalue={} {} adj=Y\n".format(astart, aend, amax, done, prefix, prefix, pval, tfarg))
    subprocess.call("chmod +x " + outpath, shell=True)
    return outpath

def write_consensus_script_OLD(directory, prefix, outfile="consensus.qsub", support=1):
    outpath = directory + "/" + outfile
    newadj = "{}.supp={}".format(prefix, support) # name of adj being generated, without .adj extension
    print "Writing {}".format(outpath)
    with open(outpath, "w") as out:
        out.write("""#!/bin/bash
#SBATCH --mem-per-cpu=90G
#SBATCH --time=60:00:00

module load dibig_tools

apple.py consensus -s {} -c {}.counts.csv -d network-data.csv {}.adj {}-*.adj*
""".format(support, prefix, newadj, prefix))
    return (outfile, newadj)

def write_consensus_script(ds, outfile="consensus.qsub", support=None, network=False):
    """Write the qsub file to generate the ADJ consensus for the specified `support'
(defaulting to 1 if None). If `network' is True, also generates the network file."""
    outpath = ds['dirname'] + "/" + outfile
    with open(outpath, "w") as out:
        out.write("""#!/bin/bash
#SBATCH --mem-per-cpu=90G
#SBATCH --time=60:00:00

module load dibig_tools

""")
        if support:
            ds['adj.supp']    = "{}.supp={}.adj".format(ds['name'], support)
            ds['counts.supp'] = "{}.counts.supp={}.csv".format(ds['name'], support)
            ds['netwrk.supp'] = "{}.netwrk.supp={}.csv".format(ds['name'], support)
            nopt = "-d " + ds['netwrk.supp'] if network else ""
            out.write("""apple.py consensus -s {} -c {} {} {} {}-*.adj*\n""".format(
                    support, ds['counts.supp'], nopt, ds['adj.supp'], ds['name']))
        else:
            ds['adj']    = ds['name'] + ".adj"
            ds['counts'] = "{}.counts.csv".format(ds['name'])
            ds['netwrk'] = "{}.netwrk.csv".format(ds['name'])
            nopt = "-d " + ds['netwrk'] if network else ""
            out.write("""apple.py consensus -s 1 -c {} {} {} {}-*.adj*\n""".format(ds['counts'], nopt, ds['adj'], ds['name']))
    return outfile

def file_columns_to_list(filename, c1, c2):
    result = []
    with open(filename, "r") as f:
        for line in f:
            if len(line) > 0 and line[0] != "#":
                line = line.rstrip("\n\r").split("\t")
                result.append([float(line[c1]), float(line[c2])])
    return result

def write_filter_script(ds, outfile="filter.qsub", mi=None, summi=None, log=None):
    directory = ds.dirname
    prefix = ds.name
    outpath = directory + "/" + outfile
    if mi:
        inadj = ds.adj_supp # "{}.adj".format(prefix)
        outadj = ds['adj.mi'] = "{}.mi.adj".format(prefix)
        flag = ""
    elif summi:
        inadj = ds['adj.mi'] # "{}.mi.adj".format(prefix)
        outadj = ds['adj.summi'] = "{}.summi.adj".format(prefix)
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

def file_columns_to_list(filename, c1, c2):
    result = []
    with open(filename, "r") as f:
        for line in f:
            if len(line) > 0 and line[0] != "#":
                line = line.rstrip("\n\r").split("\t")
                result.append([float(line[c1]), float(line[c2])])
    return result

def add_missing(data, maxv):
    # Sort data by descending first element
    data = sorted(data, key=lambda d: d[0], reverse=True)
    result = []
    m = maxv
    for i in range(len(data)):
        tg = data[i][0]
        while m > tg:
            result.append([m, 0])
            m = m - 1
        result.append(data[i])
        m = m - 1
    while m > 0:
        result.append([m, 0])
        m = m - 1
    return result

def conv_and_reverse(data):
    return sorted(data, key=lambda d: d[0], reverse=True)

def compare_mi_histograms(outfile, infile1, infile2, maxv=None):
    """Infile1 = real, infile2 = random."""
    data1 = file_columns_to_list(infile2, 0, 1)
    data2 = file_columns_to_list(infile1, 0, 1)
    # print "{} values read from {}\n{} values read from {}".format(len(data1), infile1, len(data2), infile2)
    # print data1
    # print data2
    tot1 = 0
    tot2 = 0
    maxdiff = [0, 1, 0]
    if maxv:
        data1 = add_missing(data1, maxv)
        data2 = add_missing(data2, maxv)
    else:
        data1 = conv_and_reverse(data1)
        data2 = conv_and_reverse(data2)

    with open(outfile, "w") as out:
        out.write("#Idx\tRandom\tReal\tDiff\tFPR\t% Diff\n")
        for i in range(len(data1)):
            x1 = data1[i][1]
            x2 = data2[i][1]
            tot1 += x1
            tot2 += x2
            diff = tot2-tot1
            # print "{}-{} = {} ({})".format(tot1, tot2, diff, maxdiff)
            if tot2 == 0:
                fpr = 0
            else:
                fpr = 1.0 * tot1 / tot2
            if tot1 == 0:
                pdiff = 0
            else:
                pdiff = 1.0 * diff / tot1
            out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(data1[i][0], tot1, tot2, diff, fpr, pdiff))
            # raw_input()
            if diff > maxdiff[0]:
                maxdiff[0] = diff
                maxdiff[1] = data1[i][0]
                maxdiff[2] = fpr
    return maxdiff

def merge_counts(outfile, infiles):
    data = {}
    maxn = 0
    nf = len(infiles)

    for filename in infiles:
        with open(filename, "r") as f:
            f.readline()
            for line in f:
                parsed = line.rstrip("\n\r").split("\t")
                idx = int(parsed[0])
                v = int(parsed[1])
                if idx in data:
                    data[idx] += v
                else:
                    data[idx] = v
                if idx > maxn:
                    maxn = idx
    with open(outfile, "w") as out:
        out.write("#Support\tCount\n")
        for i in range(1, maxn+1):
            if i in data:
                v = data[i]
            else:
                v = 0
            out.write("{}\t{}\n".format(i, int(round(1.0 * v / nf))))
    return maxn

### The following functions reproduce the steps performed in run-step4.sh

def get_histogram_limits(filename):
    """Read a histogram file `filename' and return the smallest and largest values in the first column."""
    hmin = 0
    hmax = 0
    with open(filename, "r") as f:
        line = f.readline().split("\t")
        hmin = float(line[0])
        for line in f:
            line = line.split("\t")
            hmax = line[0]
    hmax = float(hmax)
    return (hmin, hmax)

def do_mi_histogram(ACT, d, mirange=None, overflow=False, summi=False, final=False, log=None):
    """Write and submit a script to generate the MI or sum(MI) histogram
for dataset `d', using the specified MI range `mirange' and `overflow' flag.
Sets the 'mihist' or 'summihist' property in the dataset to the name of the
generated file. Return the ID of the submitted job, or False if the job was
not submitted because not necessary."""
    directory = d.dirname
    prefix = d.name
    infile = d.path(d.adj_supp)
    if final:
        outfile = prefix + ".final.hist.csv"
        d.finalhist = outfile
    elif summi:
        outfile = prefix + ".summi.hist.csv"
        d.summihist = outfile
    else:
        outfile = prefix + ".mi.hist.csv"
        d.mihist = outfile

    cmdline = "apple.py histogram -n 1000 -o {}".format(outfile)
    if mirange:
        cmdline += " -r {} {}".format(mirange[0], mirange[1])
    if overflow:
        cmdline += " -v"
    if final:
        cmdline += " {}/{}.summi.adj".format(directory, prefix)
    elif summi:
        cmdline += " -s"
        cmdline += " {}/{}.mi.adj".format(directory, prefix)
    else:
        cmdline += " " + infile
    if log:
        log.log("Executing: {}", cmdline)
    qsubfile = d.name + ".hist.qsub"
    with open(qsubfile, "w") as out:
        out.write("""#!/bin/bash
#SBATCH --mem-per-cpu=40G
#SBATCH --time=60:00:00

module load dibig_tools
{}
""".format(cmdline))

    if True: #ACT.missingOrStale(outfile, infile):              # *** Is this correct?
        jid = ACT.submit(qsubfile, done="hist.@.done")
        return jid
    else:
        return False

def generate_mi_histograms(act, summi=False, final=False):
    LOG = act.log
    datasets = act.datasets
    minmi = 1.0e10
    maxmi = 0.0

    LOG.log("Generating MI histograms, no range")
    nhist = 0
    for d in datasets:
        if do_mi_histogram(act, d, summi=summi, final=final, log=LOG):
            nhist += 1
    act.wait(("hist.@.done", nhist))

    for d in datasets:
        if final:
            histfile = d.finalhist
        elif summi:
            histfile = d.summihist
        else:
            histfile = d.mihist
        LOG.log("Reading histogram limits from {}", histfile)
        limits = get_histogram_limits(histfile)
        minmi = min(minmi, limits[0])
        maxmi = max(maxmi, limits[1])
    LOG.log("Overall limits: {} - {}", minmi, maxmi)
    LOG.log("Generating MI histograms, fixed range")
    nhist = 0
    for d in datasets:
        if do_mi_histogram(act, d, mirange=[minmi, maxmi], overflow=True, summi=summi, final=final, log=LOG):
            nhist += 1
    act.wait(("hist.@.done", nhist))

### Try to make filters work

class GRFilter():
    act = None
    adj_in = ""
    counts_in = ""
    adj_out = ""

    def __init__(self, act):
        self.act = act

class GRFilterSupport(GRFilter):
    """Determine optimal support from 'counts' file, generate new consensus."""
    sum = 0                     # Sum of supports
    sumfpr = 0                  # Sum of FPRs
    optsupport = 0
    fprsupport = 0

    def run(self, dry=False):
        ACT = self.act
        LOG = ACT.log

        n = 0
        dr = ACT.real
        realcounts = dr.path(dr.counts)
        for d in ACT.datasets:
            if not d.real:
                outfile = "{}.vs.{}.support.csv".format(dr.name, d.name)
                shcounts = d.path(d.counts)
                supp = compare_mi_histograms(outfile, realcounts, shcounts, maxv=ACT.nrounds+1)
                LOG.log("Comparing {} and {} into {}: maxDiff={} for support={}, FPR={}", realcounts, shcounts, outfile, supp[0], supp[1], supp[2])
                self.sum += supp[1]
                self.sumfpr += supp[2]
                n += 1

        LOG.log("sumfpr = {}", self.sumfpr)

        self.optsupport = int(round(1.0 * self.sum/n))
        self.fprsupport = int(round(1.0 * self.sumfpr/n))

        LOG.log("Re-running consensus generation with optimal support {} (fpr={})", self.optsupport, self.fprsupport)

        ## rerun consensus generation with support found in previous step

        ncons = 0
        for d in ACT.datasets:
            script = d.write_consensus_script(support=self.optsupport)
            LOG.log("Generating consensus with support={} for {}", self.optsupport, d.name)
            LOG.log("Writing .adj file: {}", d.adj_supp)
            if not dry:
                os.chdir(d.dirname)
                ACT.submit(script, done="../consensus.@.done")
                os.chdir("..")
                ncons += 1
        ACT.wait(("consensus.@.done", ncons))

class GRFilterMI(GRFilter):
    """Determine optimal MI from 'counts' file, generate new consensus."""
    sum = 0
    sumfpr = 0

    def run(self, dry=False):
        """Use compare_mi_histograms to find optimal MI threshold, filter adj files with it."""
        ACT = self.act
        LOG = ACT.log

        LOG.log("Comparing real and randomized datasets to determine optimal MI")
        generate_mi_histograms(ACT)

        self.sum = 0
        self.sumfpr = 0
        n = 0
        realmihist = ACT.real.mihist
        for d in ACT.datasets:
            if not d.real:
                outfile = "{}.vs.{}.mi.hist.csv".format(ACT.real.name, d.name)
                shmihist = d.mihist
                mi = compare_mi_histograms(outfile, realmihist, shmihist)
                LOG.log("Comparing {} and {}: maxDiff={} for mi={}, FPR={}", realmihist, shmihist, mi[0], mi[1], mi[2])
                self.sum += mi[1]
                self.sumfpr += mi[2]
                n += 1

        self.optmi = self.sum / n
        self.fprmi = self.sumfpr / n

        LOG.log("MI threshold: {} (fpr={})", self.optmi, self.fprmi)

        ## Now filter the adj files with this threshold
        nfilt = 0
        for d in ACT.datasets:
            LOG.log("Filtering dataset {} with MI={}", d.name, self.optmi)
            script = d.write_filter_script(mi=self.optmi, log=LOG)
            if not dry:
                os.chdir(d.dirname)
                ACT.submit(script, done="../filter.@.done")
                os.chdir("..")
                nfilt += 1
        ACT.wait(("filter.@.done", nfilt))

class GRFilterSumMI(GRFilter):
    """Determine optimal sum(MI) from 'counts' file, generate new consensus."""
    sum = 0
    sumfpr = 0

    def run(self, dry=False):
        """Use compare_mi_histograms to find optimal sum(MI) threshold, filter adj files with it."""

        ACT = self.act
        LOG = ACT.log

        LOG.log("Comparing real and randomized datasets to determine optimal sum(MI)")
        generate_mi_histograms(ACT, summi=True)

        self.sum = 0
        self.sumfpr = 0
        n = 0
        realmihist = ACT.real.summihist
        for d in ACT.datasets:
            if not d.real:
                outfile = "{}.vs.{}.summi.hist.csv".format(ACT.real.name, d.name)
                shmihist = d.summihist
                mi = compare_mi_histograms(outfile, realmihist, shmihist)
                LOG.log("Comparing {} and {}: maxDiff={} for sumMI={}, FPR={}", realmihist, shmihist, mi[0], mi[1], mi[2])
                self.sum += mi[1]
                self.sumfpr += mi[2]
                n += 1

        self.optsummi = self.sum / n
        self.fprsummi = self.sumfpr / n

        LOG.log("SumMI threshold: {} (fpr={})", self.optsummi, self.fprsummi)

        ## Now filter the adj files again with this threshold
        nfilt = 0
        for d in ACT.datasets:
            LOG.log("Filtering dataset {} with sum(MI)={}", d.name, self.optsummi)
            script = d.write_filter_script(summi=self.optsummi, log=LOG)
            if not dry:
                os.chdir(d.dirname)
                ACT.submit(script, done="../filter.@.done")
                os.chdir("..")
                nfilt += 1
        ACT.wait(("filter.@.done", nfilt))

class GRFilterFinal(GRFilter):

    def run(self, dry=False):

        ACT = self.act
        LOG = ACT.log

        ## Write final MI histogram
        LOG.log("Comparing real and randomized datasets to generate final MI histogram.")
        generate_mi_histograms(ACT, final=True)

        realfinalhist = ACT.real.finalhist
        finalhist = ""
        for d in ACT.datasets:
            if not d.real:
                finalhist = outfile = "{}.vs.{}.final.hist.csv".format(ACT.real.name, d.name)
                shfinalhist = d.finalhist
                mi = compare_mi_histograms(outfile, realfinalhist, shfinalhist)

        ## *** TODO: average final histograms into single one, or get apple.py to average them.

        ## Write stats on all adj files after each filtering step
        all_adj_stats(ACT, "adj-stats.csv")
        with open("adj-stats.csv", "r") as f:
            stats = f.read()
        LOG.log("ADJ stats after each filtering step:\n{}", stats)

        ## Convert final adj file to cytoscape
        final = ACT.real.path(ACT.real.adj_summi)
        cyto = ACT.real.name + ".cyto.csv"
        LOG.log("Converting final file {} to cytoscape file {}", final, cyto)

        ## *** Hack: we're simply taking a random finalhist. We should use the averaged one from above.
        ACT.shell("module load dibig_tools; apple.py convert ac -f {} {} {}".format(finalhist, final, cyto))
        ## Sort the resulting cytoscape file, otherwise later conversion to ADJ will fail.
        ACT.shell("rm -f unsorted; mv {} unsorted; head -1 unsorted > {}; tail -n +2 unsorted | sort >> {}".format(cyto, cyto, cyto))

        ## If translation file specified, do translation

        if ACT.translation:
            cytopre = cyto
            cyto = ACT.real.name + ".tr.csv"
            LOG.log("Translating gene names in {} to {} using map ../{}", cytopre, cyto, ACT.translation)
            ACT.shell("module load dibig_tools; apple.py translate ../{} {} {}".format(ACT.translation, cytopre, cyto))
            finalpre = final
            final = ACT.real.name + ".tr.adj"
            LOG.log("Converting translated file {} to {}", cyto, final)
            ACT.shell("module load dibig_tools; apple.py convert ca {} {}".format(cyto, final))

        ## Produce connections file and CX file, and print top genes.
        connections = ACT.real.name + ".conn.csv"
        ACT.shell("module load dibig_tools; apple.py convert co {} {}".format(cyto, connections))

        cxfile = ACT.real.name + ".cx"
        afile = CXattributes(ACT)
        dbfile = "-n " + ACT.genesdb if ACT.genesdb else ""
        ACT.shell("module load dibig_tools; apple.py convert cx {} {} {} {}".format(afile, dbfile, cyto, cxfile))
        ACT.shell("cut -f 1,2 {} | head -{} > tophubs.txt".format(connections, ACT.tophubs))
        with open("tophubs.txt", "r") as f:
            tophubs = f.read()

        LOG.log("Top {} hubs by number of connections:\n".format(ACT.tophubs) + tophubs)

        LOG.log("Writing tophubnet.cy")
        ACT.shell("module load dibig_tools; apple.py extract -a -c -o tophubnet.cy {} tophubs.txt".format(cyto))
        cxfile = ACT.real.name + "-tophubs.cx"
        LOG.log("Converting tophubnet.cy to CX format file {}", cxfile)

        ACT.shell("module load dibig_tools; apple.py convert cx {} {} tophubnet.cy {}".format(afile, dbfile, cxfile))

        return True

def CXattributes(ACT):
    if len(ACT.attributes) > 0:
        with open("attributes.txt", "w") as out:
            for a in ACT.attributes:
                out.write("{}: {}\n".format(a[0], a[1]))
            out.write("version: {}\n".format(date.today().isoformat()))
            out.write("networkType: Genetic interactions\n")
        return "-a attributes.txt"
    else:
        return ""

### Write adj stats report

def all_adj_stats(ACT, outfile):
    cmdline = "module load dibig_tools; apple.py stats -o {} ".format(outfile)
    for d in ACT.datasets:
        cmdline += " " + d.dirname + "/" + d.adj_supp
    for d in ACT.datasets:
        cmdline += " " + d.dirname + "/" + d.adj_mi
    for d in ACT.datasets:
        cmdline += " " + d.dirname + "/" + d.adj_summi
    ACT.shell(cmdline)

### Temporary

def swap_columns(infile, outfile, col1, col2):
    with open(outfile, "w") as out:
        with open(infile, "r") as f:
            out.write(f.readline())
            for line in f:
                parsed = line.split("\t")
                x = parsed[col1]
                parsed[col1] = parsed[col2]
                parsed[col2] = x
                out.write("\t".join(parsed))

### Clean input expressions file

def removeZeros(infile, outfile, prop=0.5, genecols=2):
    """Remove lines from `infile' in which the proportion of zeros is equal to or higher than `prop'. `genecols' is the number of columns containing gene identifiers at the beginning of each row. Writes filtered lines to `outfile'."""
    nin = 0
    nout = 0
    with open(infile, "r") as f:
        hdr = f.readline()
        columns = hdr.split("\t")
        ncols = len(columns)-genecols
        maxzeros = ncols*prop
        with open(outfile, "w") as out:
            out.write(hdr)
            while True:
                line = f.readline()
                if line == '':
                    break
                nin += 1
                pline = line.rstrip("\r\n").split("\t")
                nzeros = 0
                for v in pline[genecols:]:
                    if float(v) == 0:
                        nzeros += 1
                if nzeros < maxzeros:
                    out.write(line)
                    nout += 1
    return (nin, nout)
