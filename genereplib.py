import os
import glob
import os.path
import subprocess

### Functions implementing script steps

def submit_aracne_array(infile, outfile="aracne-all.sh", tfs=False, astart=0, aend=99, amax=40, pval="1e-7", dpi="0.1", done=None):
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
        out.write("submit -t {}-{}%{} {} aracne-array.qsub {} {} pval={} dpi={}".format(astart, aend, amax, done, prefix, prefix, pval, dpi))
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

def write_consensus_script(directory, prefix, outfile="consensus.qsub", support=1):
    outpath = directory + "/" + outfile
    print "Writing {}".format(outpath)
    with open(outpath, "w") as out:
        out.write("""#!/bin/bash
#SBATCH --mem-per-cpu=80G
#SBATCH --time=60:00:00

module load dibig_tools

apple.py consensus -s {} -c {}.counts.csv -d network-data.csv {}.adj {}-*.adj*
""".format(support, prefix, prefix, prefix))
    return outfile

def file_columns_to_list(filename, c1, c2):
    result = []
    with open(filename, "r") as f:
        for line in f:
            if len(line) > 0 and line[0] != "#":
                line = line.rstrip("\n\r").split("\t")
                result.append([float(line[c1]), float(line[c2])])
    return result

def write_filter_script(directory, prefix, outfile="filter.qsub", mi=None, summi=None):
    outpath = directory + "/" + outfile
    if mi:
        inadj = "{}.adj".format(prefix)
        outadj = "{}.mi.adj".format(prefix)
        flag = ""
    elif summi:
        inadj = "{}.mi.adj".format(prefix)
        outadj = "{}.summi.adj".format(prefix)
        flag = "-t"
        mi = summi
    else:
        print "Error: one of `mi' and `summi' must be specified!"
        sys.exit(1)
    print "Writing {}".format(outpath)
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
            if tot2 == 0:
                fpr = 0
            else:
                fpr = 1.0 * tot1 / tot2
            if tot1 == 0:
                pdiff = 0
            else:
                pdiff = 1.0 * diff / tot1
            out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(data1[i][0], tot1, tot2, diff, fpr, pdiff))
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

def do_mi_histogram(act, d, mirange=None, overflow=False, summi=False):
    """Write and submit a script to generate the MI or sum(MI) histogram
for dataset `d', using the speficied MI range `mirange' and `overflow' flag.
Sets the 'mihist' or 'summihist' property in the dataset to the name of the
generated file. Return the ID of the submitted job."""
    directory = d['dirname']
    prefix = d['name']
    if summi:
        outfile = prefix + ".summi.hist.csv"
        d['summihist'] = outfile
    else:
        outfile = prefix + ".mi.hist.csv"
        d['mihist'] = outfile
    cmdline = "apple.py histogram -n 1000 -o {}".format(outfile)
    if mirange:
        cmdline += " -r {} {}".format(mirange[0], mirange[1])
    if overflow:
        cmdline += " -v"
    if summi:
        cmdline += " -s"
        cmdline += " {}/{}.mi.adj".format(directory, prefix)
    else:
        cmdline += " {}/{}.adj".format(directory, prefix)
    qsubfile = prefix + ".hist.qsub"
    with open(qsubfile, "w") as out:
        out.write("""#!/bin/bash
#SBATCH --mem-per-cpu=40G
#SBATCH --time=60:00:00

module load dibig_tools
{}
""".format(cmdline))
    jid = act.submit(qsubfile, done="hist.@.done")
    return jid

def generate_mi_histograms(act, summi=False):
    LOG = act.log
    datasets = act.datasets
    minmi = 1.0e10
    maxmi = 0.0

    LOG.log("Generating MI histograms, no range")
    nhist = 0
    for d in datasets:
        do_mi_histogram(act, d, summi=summi)
        nhist += 1
    act.wait(("hist.@.done", nhist))

    for d in datasets:
        if summi:
            histfile = d['summihist']
        else:
            histfile = d['mihist']
        LOG.log("Reading histogram limits from {}", histfile)
        limits = get_histogram_limits(histfile)
        minmi = min(minmi, limits[0])
        maxmi = max(maxmi, limits[1])
    LOG.log("Overall limits: {} - {}", minmi, maxmi)
    LOG.log("Generating MI histograms, fixed range")
    nhist = 0
    for d in datasets:
        do_mi_histogram(act, d, mirange=[minmi, maxmi], overflow=True, summi=summi)
        nhist += 1
    act.wait(("hist.@.done", nhist))

### Write adj stats report

def all_adj_stats(ACT, outfile):
    cmdline = "module load dibig_tools; apple.py stats -o {} ".format(outfile)
    for d in ACT.datasets:
        cmdline += " " + d['dirname'] + "/" + d['name'] + ".adj"
    for d in ACT.datasets:
        cmdline += " " + d['dirname'] + "/" + d['name'] + ".mi.adj"
    for d in ACT.datasets:
        cmdline += " " + d['dirname'] + "/" + d['name'] + ".summi.adj"
    print "Executing: " + cmdline
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
