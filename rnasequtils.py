# (c) 2016, A. Riva, DiBiG
# RNA-Seq Analysis pipeline helper functions, v2.0

import math
import Utils
import Plots

def writeMatrixScript(scriptfile, samples, outfile, mode='g', erccdb=False, column=None):
    """Write a qsub script to call 'rnaseqtools.py matrix' on the specified samples. Reads RSEM *.genes.results files
if mode is 'g', RSEM *.isoforms.results files if mode is 'i', and kallisto abundance.tsv files if mode is 'k'. 
Extracts column `expected_count' (or 'tpm' for kallisto) by default unless a different column is specified.
If `erccdb' is specified it should be a pathname to and ERCC databases. In this case, ERCC normalization is
performed."""
    with open(scriptfile, "w") as out:
        out.write("""#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=10G

module load rsem/1.2.29
module load dibig_tools

rnaseqtools.py matrix """)
        if column:
            out.write(" -c " + column)
        elif mode in ['g', 'i']:
            out.write(" -c expected_count")
        elif mode == 'k':
            out.write(" -c tpm")
        if erccdb:
            out.write(" -e -ercc " + erccdb + " -mix " + ",".join([s['mix'] for s in samples]))
        for s in samples:
            if mode == 'g':
                out.write(" {}.genes.results".format(s['name']))
            elif mode == 'i':
                out.write(" {}.isoforms.results".format(s['name']))
            elif mode == 'k':
                out.write(" {}.kallisto.d/abundance.tsv".format(s['name']))
        out.write(" > {}\n\n".format(outfile))
    return scriptfile

def writeContrastMatrix(contrast, samples, ntest, nctrl, fdr, erccdb=None, mode="g"):
    with open(contrast[mode+'rsemqsub'], "w") as out:
        out.write("""#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=10G

module load rsem/1.2.29
module load dibig_tools

#rsem-generate-data-matrix 

rnaseqtools.py matrix """)
        if erccdb:
            out.write(" -e -ercc " + erccdb + " -mix " + ",".join([s['mix'] for s in samples]))
        for s in samples:
            if mode == 'g':
                out.write(" {}.genes.results".format(s['name']))
            else:
                out.write(" {}.isoforms.results".format(s['name']))
        out.write(" > {}\n\n".format(contrast[mode+'matrix']))
        out.write("rsem-run-ebseq {} {},{} {}\n".format(contrast[mode+'matrix'], ntest, nctrl, contrast[mode+'diff']))
        # out.write("rsem-control-fdr --soft-threshold {} {} {}\n".format(contrast[mode+'diff'], fdr, contrast[mode+'fdr']))  # Soft threshold is confusing...
        out.write("rsem-control-fdr {} {} {}\n".format(contrast[mode+'diff'], fdr, contrast[mode+'fdr']))
    return contrast[mode+'rsemqsub']

def filterDiff(infile, outfile, fc, translation=None, wanted=[], wantedNames=[], biotype=None):
    """Read RSEM differential expression results from `infile' and write them to `outfile', converting
the fold changes to log2 and discarding records having log2(FC) < `fc'. If `translation' is provided
it should be a dictionary mapping gene or transcript IDs (first column) to additional annotations. The
`extra' list provides the header names for the additional annotations. If `biotype' is specified, the
annotations are assumed to contain a 'biotype' as the second element, and the output file will contain only records having
the specified biotype. Returns a tuple containing the number of records written, number of overexpressed,
number of underexpressed, and a dictionary mapping each biotype to is number of occurrences."""
    data = []
    biotypes = {}
    missing = ["???"]*len(wanted)
    ng = nup = ndown = 0
    with open(infile, "r") as f:
        f.readline()
        for line in f:
            good = True
            parsed = line.rstrip("\r\n").split("\t")
            log2fc = math.log(float(parsed[3]), 2)
            if abs(log2fc) >= fc:
                if translation:
                    k = parsed[0].strip('"')
                    if k in translation:
                        trans = translation[k]
                        if biotype:
                            biot = trans['gene_biotype']
                            if biotype != biot:
                                good = False
                            else:
                                Utils.dinc(biot, biotypes)
                    else:
                        good = False
                    if good:
                        toAdd = [ trans[tag] for tag in wanted ]
                    else:
                        toAdd = missing
                    row = [log2fc] + toAdd + [str(log2fc), parsed[1], parsed[5], parsed[6]]
                else:
                    row = [log2fc, parsed[0], str(log2fc), parsed[1], parsed[5], parsed[6]]

                if good:
                    ng += 1
                    if log2fc > 0:
                        nup += 1
                    else:
                        ndown += 1
                    data.append(row)

    data.sort(key=lambda d: d[0], reverse=True)

    with open(outfile, "w") as out:
        if wanted:
            out.write("\t".join(wantedNames) + "\tLog2(FC)\tP-value\tExp0\tExp1\n")
        else:
            out.write("Gene\tLog2(FC)\tP-value\tExp0\tExp1\n")
        for d in data:
            out.write("\t".join(d[1:]) + "\n")
    return (ng, nup, ndown, biotypes)

### Scatterplot for two columns of expression values

def avgCols(row, columns, nc):
    v = 0
    for c in columns:
        v += float(row[c])
    return v/nc

def exprScatterplot(filename, exprfile, cols1, cols2, label1=None, label2=None, fc=1):
    data = []
    nc1 = len(cols1)
    nc2 = len(cols2)
    with open(exprfile, "r") as f:
        hdr = f.readline().rstrip("\r\n").split("\t")
        label1 = label1 or hdr[cols1[0]]
        label2 = label2 or hdr[cols2[0]]
        for line in f:
            parsed = line.rstrip("\r\n").split("\t")
            v1 = avgCols(parsed, cols1, nc1)
            v2 = avgCols(parsed, cols2, nc2)
            data.append([ v1, v2 ])

    p = Plots.Scatterplot(data=data, title='Scatterplot', ylabel=label1, xlabel=label2)
    p.fc = fc
    p.plot(filename)
    return p
