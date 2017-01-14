# (c) 2015, A. Riva, DiBiG
# Differential Methylation Analysis pipeline, v2.0

from glob import glob

import Utils

def initializeContrasts(ACT):
    for c in ACT.sc.contrasts:
        label = c['label']
        c['diff'] = ACT.getConf("diff", section=label, default=ACT.diff)
        c['pval'] = ACT.getConf("pval", section=label, default=ACT.pval)
        c['genesdb'] = ACT.getConf("genesdb", section=label, default=ACT.genesdb)
        c['genes'] = ACT.getConf("genes", section=label, default=ACT.genes)
        c['regions'] = ACT.getConf("regions", section=label, default=ACT.regions)
        c['size'] = ACT.getConf("size", section=label, default=ACT.size)
        c['mode'] = ACT.getConf("mode", section=label, default=ACT.mode)
        c['nsites'] = ACT.getConf("nsites", section=label, default=ACT.nsites)

        ## Fix and validate values
        if c['diff'] != None:
            c['diff'] = float(c['diff'])
        if c['pval'] != None:
            c['pval'] = float(c['pval'])
        if c['genesdb'] != None:
            c['genesdb'] = ACT.checkPath(c['genesdb'])

        if c['genes'] == None:
            c['genes'] = []
        else:
            c['genes'] = c['genes'].split(",")

        if c['regions'] != None:
            for r in c['regions']:
                if not r in "bipedBIPED":
                    print "Error: `regions' should only contain P, B, E, I, or D."
                    sys.exit()
                    
        if c['size'] != None:
            c['size'] = int(c['size'])

        if c['mode'] != None:
            if not c['mode'] in ['avg', 'abs', 'max', 'min', 'bal']:
                print "Error: `mode' should be one of avg, abs, max, min, or bal."
                sys.exit()
        
        if c['nsites'] != None:
            c['nsites'] = int(c['nsites'])

def writeBSMAPscripts(ACT):
    """Write a qsub script for each sample. The script will align each
fastq or pair of fastqs with bsmap. At the end, all resulting BAM files
are merged."""
    bsmap_bin = "/apps/dibig_tools/dmap/bsmap"

    for smp in ACT.sc.samples:
        print "writing bsmap for {}".format(smp['name'])
        smp['bsmapscript'] = "bsmap-{}.qsub".format(smp['name'])
        smp['bam'] = smp['name'] + ".bam"
        ACT.log.log("Writing BSMAP script {}", smp['bsmapscript'])
        rrbsopt = ""
        if ACT.rrbs:
            rrbsopt = "-D " + ACT.rrbs
        with open(smp['bsmapscript'], "w") as out:
            out.write("""#!/bin/bash

#SBATCH --mem=20G
#SBATCH --time=60:00:00
#SBATCH --ntasks=4

module load samtools
module load bamtools

""")
            for rs in smp['readsets']:
                if rs['paired']:
                    out.write("{} -d {} -a {} -b {} -o {} -p 4 {}\n".format(bsmap_bin, ACT.reference, rs['left'], rs['right'], rs['name'] + ".bam", rrbsopt))
                else:
                    out.write("{} -d {} -a {} -o {} -p 4 {}\n".format(bsmap_bin, ACT.reference, rs['left'], rs['name'] + ".bam", rrbsopt))

            if len(smp['readsets']) == 1:
                out.write("mv {}.bam {}.bam\n".format(smp['readsets'][0]['name'], smp['name']))
            else:
                out.write("bamtools merge -out {}.bam".format(smp['name']))
                for rs in smp['readsets']:
                    out.write(" -in {}.bam".format(rs['name']))
                out.write("\n")
                
def writeMMAPconfs(ACT):
    bsmap_bin = "/apps/dibig_tools/dmap/bsmap"

    for smp in ACT.sc.samples:
        smp['mmapconf'] = confname = "mmap-{}.conf".format(smp['name'])

        ACT.log.log("Writing MMAP configuration file {}", confname)
        with open(confname, "w") as out:
            out.write("""MOABS_BIN=$HPC_MOABS_BIN

[INPUT]
""")
            for rs in smp['readsets']:
                out.write("{}_1={}\n".format(rs['name'], rs['left']))
                if rs['paired']:
                    out.write("{}_2={}\n".format(rs['name'], rs['right']))
            out.write("""
[TASK]
Program=MMAP
Label={}
Parallel=THREAD

[MMAP]
Path={}
d={}
""".format(smp['name'], bsmap_bin, ACT.reference))
            if ACT.parallel:
                out.write("p={}\n".format(ACT.parallel))

def writeMCOMPqsub(ACT):
    qsubnames = []
    idx = 1
    directives = "#SBATCH --mem=20G\n#SBATCH --time=60:00:00\n#SBATCH --nodes=1\n#SBATCH --ntasks={}\n".format(ACT.parallel)

    for contr in ACT.sc.contrasts:
        test = contr['test']
        ctrl = contr['control']
        testcond = ACT.sc.findCondition(test)
        controlcond = ACT.sc.findCondition(ctrl)
        contr['mcompout'] = "{}.vs.{}.mcomp.csv".format(test, ctrl)
        qsub = "mcomp-{}.qsub".format(idx)
        idx += 1
        qsubnames.append(qsub)
        print "Writing MCOMP submission script {}".format(qsub)
        with open(qsub, "w") as out:
            out.write("""#!/bin/bash
{}
module load moabs
mcomp -r {} -r {} -c {} --doDmcScan 0 --doDmrScan 0 -p {}
""".format(directives, controlcond['callbed'], testcond['callbed'], contr['mcompout'], ACT.parallel))
    return qsubnames

### Differential methylation analysis

# TAGS = ['insigOrLowDepth', 'StrongHypo', 'Hypo', 'Hyper', 'StrongHyper']
TAGS = ['.', '--', '-', '+', '++']

def callDmeth(ACT, thisDiff, thisPval):
    if thisPval <= ACT.pval:
        if thisDiff > 0:
            if thisDiff > ACT.diff:
                # return "StrongHyper"
                return '++'
            else:
                # return "Hyper"
                return '+'
        else:
            if thisDiff < -ACT.diff:
                # return "StrongHypo"
                return '--'
            else:
                # return "Hypo"
                return '-'
    else:
        # return "insigOrLowDepth"
        return '.'

def diffmethSummary(ACT, prefix, dry):
    result = []
    notsig = TAGS[0]
    LOG = ACT.log

    nc = 0

    for contr in ACT.sc.contrasts:
        mc = contr['mcompout']
        ctrl = contr['control']
        test = contr['test']
        fullfile = "{}.vs.{}.full.csv".format(test, ctrl)
        sigfile = "{}.vs.{}.sig.csv".format(test, ctrl)
        fullfilex = "{}.vs.{}.full.xlsx".format(test, ctrl)
        sigfilex = "{}.vs.{}.sig.xlsx".format(test, ctrl)
        bedfile = "{}.vs.{}.bed".format(test, ctrl)
        data = {'test': test, 'control': ctrl, 'total': 0, 'filename': mc, 
                'fullfile': fullfile, 'sigfile': sigfile, 'fullfilex': fullfilex, 'sigfilex': sigfilex, 'bedfile': bedfile}
        for t in TAGS:
            data[t] = 0
        with open(fullfile, "w") as full:
            with open(sigfile, "w") as sig:
                full.write("Chrom\tPos\t{}\t{}\tDiff\tPval\tCall\n".format(ctrl, test))
                sig.write("Chrom\tPos\t{}\t{}\tDiff\tPval\tCall\n".format(ctrl, test))
                LOG.log("Analyzing {} with pval={}, diff={}", mc, ACT.pval, ACT.diff)
                nin = 0
                nout = 0
                with open(mc, "r") as f:
                    f.readline()        # skip header
                    for line in f:
                        line = line.rstrip("\n").split("\t")
                        if line[12] != "NA":
                            nin += 1
                            data['total'] += 1
                            thisDiff = float(line[12])
                            thisPval = float(line[18])
                            call = callDmeth(ACT, thisDiff, thisPval)
                            data[call] += 1
                            outl = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line[0], line[1], line[5], line[9], thisDiff, thisPval, call)
                            full.write(outl)
                            if call != notsig:
                                sig.write(outl)
                                nout += 1
        LOG.log("{}: {} sites, {} differentially methylated.", mc, nin, nout)
        title = "{}.vs.{}-dmc".format(test, ctrl)
        desc  = "Differential methylation in {} vs. {}".format(test, ctrl)

        qsub = "{}-{}.qsub".format(prefix, nc)
        with open(qsub, "w") as out:
            out.write("""#!/bin/bash

#SBATCH --mem=5G
#SBATCH --time=1:00:00

module load dibig_tools
csvtoxls.py {} -q {}
csvtoxls.py {} -q {}
bamToWig.py -m -o {} -t {} -d "{}" {}
""".format(data['fullfilex'], data['fullfile'], data['sigfilex'], data['sigfile'], data['bedfile'], title, desc, data['sigfile']))
            
        if not dry:
            ACT.submit(qsub, done="diffmeth.@.done")
            nc += 1
        result.append(data)
    ACT.wait(("diffmeth.@.done", nc))

    LOG.log("Deleting temporary diffmeth files.")
    ACT.shell("rm -f {}*".format(prefix))
    return result

def getCscallReport(cond):
    """Returns the 'Total' line from the ...-report.csv file for condition `cond'."""
    with open(cond['callreport'], "r") as f:
        for line in f:
            if line.startswith("Total:"):
                parsed = line.rstrip("\r\n").split("\t")
                return parsed[1:]
    return None

### Plotting

def drawBEDplots(ACT, cond, dry):
    ACT.mkdir("plots")
    name = cond['name']
    htmlfile = "{}-metyhlation.html".format(name)
    pngfilelist = "pngfiles-" + name
    cond['methplots'] = htmlfile
    bedfile = cond['callbed']
    outpatt = "plots/methprofile-" + name + "-{}.png"
    title = "Methylation profile - {}"
    if not dry:
        ACT.shell('module load dibig_tools; plotter -o {} -f $HPC_DIBIGTOOLS_DIR/lib/fonts/liberation/ bed -i {} -title "{}" -xlabel "Chromosome position" -ylabel "% methylation" > {}', 
                  outpatt, bedfile, title, pngfilelist)

    with open(pngfilelist, "r") as f:
        plotfiles = f.readlines()
    plotfiles = [ s.strip("\n").split("\t") for s in plotfiles ]

    with open(htmlfile, "w") as out:
        out.write("""<!DOCTYPE html>
<HTML>
<HEAD>
<TITLE>Methylation level plots - {}</TITLE>
<STYLE>
BODY {{
  background: #ffff;
}}
IMG {{
  border: 2px solid grey;
}}
</STYLE>
</HEAD>
<BODY>
<CENTER>
""".format(name))
        out.write("<h1>Sample {}</h1><br>\n".format(name))
        for p in plotfiles:
            # out.write("<H3>{}</H3>".format(p[0]))
            out.write("<IMG src='{}'/><BR><BR>\n".format(p[1]))
    return htmlfile

def drawScatterplots(ACT):
    htmlfile = "scatterplots.html"
    with open(htmlfile, "w") as out:
        out.write("<!DOCTYPE html><HTML><HEAD><TITLE>Methylation scatterplots</TITLE></HEAD>\n<BODY bgcolor=#ffffff>\n<CENTER>\n")
        for contr in ACT.sc.contrasts:
            mc = contr['mcompout']
            ctrl = contr['control']
            test = contr['test']
            png = "plots/{}.vs.{}.scatt.png".format(test, ctrl)
            contr['scatterplot'] = png
            ACT.shell("module load dibig_tools; plotter -o {} dscatt -i {} -width {} -height {} -title 'Methylation level - {} vs {}' -xlabel {} -ylabel {}".format(
                    png, mc, ACT.scatterplot_size, ACT.scatterplot_size, test, ctrl, ctrl, test))
            out.write("<IMG src='{}'/><BR><BR>\n".format(png))
        out.write("</CENTER>\n</BODY>\n</HTML>\n")
    return htmlfile

def drawDMCplots(ACT, dry):
    for d in ACT.diffmeth:
        ctrl = d['control']
        test = d['test']
        sig = d['sigfile']
        png = "plots/" + test + ".vs." + ctrl + "-{}.dmc.png" # png filename with placeholder for chrom        
        htmlfile = "{}.vs.{}.dmc.html".format(test, ctrl)
        d['dmcplot'] = png
        d['dmcfile'] = htmlfile
        title = "Differential methylation - {}"
        if not dry:
            ACT.shell('module load dibig_tools; plotter -o {} -f $HPC_DIBIGTOOLS_DIR/lib/fonts/liberation/ dmc -i {} -xlabel "Chromosome position" -ylabel "% methylation" -title "{}" > pngfiles2',
                      png, sig, title)

        with open("pngfiles2", "r") as f:
            plotfiles = f.readlines()
        plotfiles = [ s.strip("\n").split("\t") for s in plotfiles ]

        with open(htmlfile, "w") as out:
            out.write("<!DOCTYPE html><HTML><HEAD><TITLE>Differential methylation - {} vs {}</TITLE></HEAD>\n<BODY bgcolor=#ffffff>\n<CENTER>\n".format(test, ctrl))
            out.write("<h1>Differential methylation - {} vs {}</h1><br>\n".format(test, ctrl))
            for p in plotfiles:
                out.write("<IMG src='{}'/><BR><BR>\n".format(p[1]))

def runAvgMeth(ACT, cond1, cond2):
    """Run the dmaptools avgmeth command con the raw files for the two conditions `cond1' and `cond2', returning a dictionary with the results."""
    result = {}

    SC = ACT.sc
    cnd1 = SC.findCondition(cond1)
    cnd2 = SC.findCondition(cond2)
    out = "avgmeth.txt"
    cmdline = "module load dibig_tools; dmaptools.py avgmeth {} {} {}".format(cnd1['matreport'], cnd2['matreport'], out)
    ACT.shell(cmdline)
    if ACT.checkFile(out):
        with open(out, "r") as f:
            for line in f:
                line = line.rstrip("\r\n")
                p = line.find(":")
                if p > 0:
                    result[line[:p]] = line[p+1:]
    return result

def runMethHist(ACT, cond1, cond2, outfile):
    SC = ACT.sc
    cnd1 = SC.findCondition(cond1)
    cnd2 = SC.findCondition(cond2)
    cmdline = "module load dibig_tools; dmaptools.py histmeth {} {} {}".format(cnd1['matreport'], cnd2['matreport'], outfile)
    ACT.shell(cmdline)
    if ACT.checkFile(outfile):
        return Utils.fileToList(outfile)
    else:
        return None
