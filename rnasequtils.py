# (c) 2016, A. Riva, DiBiG
# RNA-Seq Analysis pipeline helper functions, v2.0

import math
import Utils

def writeMatrixScript(scriptfile, samples, outfile, mode='g', erccdb=False, column=None):
    """Write a qsub script to call 'rnaseqtools.py matrix' on the specified samples. Gene level if mode is 'g',
isoform level if mode is 'i'. Extracts column `expected_count' by default unless a different column is specified.
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
        if erccdb:
            out.write(" -e -ercc " + erccdb + " -mix " + ",".join([s['mix'] for s in samples]))
        for s in samples:
            if mode == 'g':
                out.write(" {}.genes.results".format(s['name']))
            elif mode == 'i':
                out.write(" {}.isoforms.results".format(s['name']))
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
                        biot = trans['gene_biotype']
                        if biotype:
                            if biotype != biot:
                                good = False
                            else:
                                Utils.dinc(biot, biotypes)
                    else:
                        trans = missing
                    if good:
                        toAdd = [ trans[name] for name in wanted ]
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

### Generation of hub for RNA-seq data

