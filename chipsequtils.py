# (c) 2017, A. Riva, DiBiG
# ChIP-Seq Analysis pipeline helper functions, v2.0

import Utils
import Plots

def readHomerAnnots(filename):
    (data, header) = Utils.fileToList(filename, hdr=True)
    data2 = []
    for row in data:
        if row[0] == "Annotation":
            break
        data2.append(row)
    data2.sort(key=lambda x: float(x[1]), reverse=True)
    data2.insert(0, header)
    return data2

def writeClassifyScript(conditions, genesdb):
    """Write a qsub script to classify all peaks in all conditions."""
    scriptfile = "classify.qsub"
    with open(scriptfile, "w") as out:
        out.write("""#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=10G

module load dibig_tools

DB={}

""".format(genesdb))
        for c in conditions:
            infile = c['tags'] + "/peaks.txt"
            outfile = c['peaksdistr'] = c['name'] + "-peaks-distr.xlsx"
            out.write("files=$(genes.py split -x {} -db $DB {}:2)\n".format(c['name'], infile))
            out.write("csvtoxls.py {} $files\n\n".format(outfile))

            infile = c['tags'] + "/regions.txt"            
            outfile = c['regionsdistr'] = c['name'] + "-regions-distr.xlsx"
            out.write("files=$(genes.py split -x {} -db $DB {}:2)\n".format(c['name'], infile))
            out.write("csvtoxls.py {} $files\n\n".format(outfile))

            infile = c['tags'] + "/superEnhancers.txt"
            outfile = c['enhancersdistr'] = c['name'] + "-enhancers-distr.xlsx"
            out.write("files=$(genes.py split -x {} -db $DB {}:2)\n".format(c['name'], infile))
            out.write("csvtoxls.py {} $files\n\n".format(outfile))
    return scriptfile

def classificationBarChart(filename, conditions, what, title):
    """Produce a bar chart showing the classification for the sites of type `what' (either
peaks, regions, or superEnhancers)."""
    groups = ['Upstream', 'Exon', 'CodingExon', 'Intron', 'Downstream', 'Intergenic']
    names = []
    data = []
    for c in conditions:
        names.append(c['name'])
        infile = c['tags'] + "/" + what + "-Summary.csv"
        values = Utils.fileToDict(infile, toInt=False, column=2)
        row = [ values[g] for g in groups ]
        row = [ float(x[:-1]) for x in row]
        data.append(row)
    p = Plots.BarChart(data=data, series=names, xticklabels=groups, ylabel="Percentage", title=title)
    p.plot(filename)

