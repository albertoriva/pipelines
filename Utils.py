# (c) 2016, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import os
import os.path
import csv
import gzip
import time
import pysam
import string
import random
import subprocess

def genOpen(filename, mode):
    """Generalized open() function - works on both regular files and .gz files."""
    (name, ext) = os.path.splitext(filename)
    if ext == ".gz":
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def dget(key, dictionary, default=None):
    """Return the value associated with `key' in `dictionary', if present, or `default'."""
    if key in dictionary:
        return dictionary[key]
    else:
        return default

def dinc(key, dictionary, default=0):
    """Increment the value associated with `key' in dictionary. If key is not present, initialize it with `default'."""
    if key in dictionary:
        dictionary[key] += 1
    else:
        dictionary[key] = default + 1
    return dictionary[key]

def id_generator(size=6, chars=string.ascii_uppercase + string.digits, prefix=''):
    """From: http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python"""
    return prefix + ''.join(random.choice(chars) for _ in range(size))

def safeInt(v, default=None):
    try:
        return int(v)
    except:
        return default

def safeReadIntFromFile(filename, default=None, maxtries=20, delay=1):
    """Read an integer from the first line of file `filename', using `default' if the
line does not contain an integer. If the line is empty, wait `delay' seconds and
retry, up to `maxtries' times."""
    while True:
        with open(filename, "r") as f:
            c = f.readline().rstrip("\n")
        if c == '':
            maxtries -= 1
            if maxtries == 0:
                return default
            else:
                time.sleep(delay)
        else:
            return safeInt(c, default)

def safeReadLineFromFile(filename, default=None, maxtries=20, delay=1):
    """Read and return the first line of file `filename'. If the line is empty, wait `delay' seconds and
retry, up to `maxtries' times, after which return `default'."""
    while True:
        with open(filename, "r") as f:
            c = f.readline().rstrip("\n")
        if c == '':
            maxtries -= 1
            if maxtries == 0:
                return default
            else:
                time.sleep(delay)
        else:
            return c

def safeRemoveFile(filename):
    if os.path.isfile(filename):
        try:
            os.remove(filename)
        except OSError:
            pass

def filenameNoExt(s):
    return os.path.splitext(os.path.basename(s))[0]

def parseLine(line):
    return line.rstrip("\r\n").split("\t")

def fmt(n):
    return "{:,}".format(n)

def pct(x, den):
    return "{:.2f}%".format(x*100.0/den)

def f2dd(x):
    return "{:.2f}".format(x)

def f3dd(x):
    return "{:.3f}".format(x)

def up(x):
    return "<span class='upreg'>{}</span>".format(x)

def UP(x):
    return "<span class='upreg'><b>{}</b></span>".format(x)

def down(x):
    return "<span class='dnreg'>{}</span>".format(x)

def DOWN(x):
    return "<span class='dnreg'><b>{}</b></span>".format(x)

def plural(x):
    return "" if x == 1 else "s"

def fileToDict(filename, toInt=True, column=1):
    """Read `filename' and return a dictionary having as keys the strings
in the first column and as values the contents of the specified `column' 
(converted to int if `toInt' is True)."""
    result = {}
    with open(filename, "r") as f:
        for line in f:
            parsed = line.rstrip("\r\n").split("\t")
            if toInt:
                result[parsed[0]] = int(parsed[column])
            else:
                result[parsed[0]] = parsed[column]
    return result

def fileToList(filename, delimiter='\t', hdr=True):
    """Read delimited file `filename' and return its contents as a list of lists. If `hdr' is True (the
default) treat the first line as header. Returns a tuple (data, header)."""
    result = []
    header = []
    with open(filename, "r") as f:
        r = csv.reader(f, delimiter=delimiter)
        if hdr:
            header = r.next()
        for line in r:
            result.append(line)
    return (result, header)

def linkify(url, text=None, target=None):
    if not text:
        text = url
    if target:
        tg = " target='{}'".format(target)
    else:
        tg = ""
    return "<A href='{}'{}>{}</A>".format(url, tg, text)

### Should this be here?

def detectFileFormat(filename):
    with genOpen(filename, "r") as f:
        line = f.readline()
        if len(line) > 0:
            if line[0] == ">":
                return "fasta"
            elif line[0] == "@":
                return "fastq"
        return "?"

def fastqcPath(basedir, f):
    base = os.path.split(f)[1]
    if base.endswith(".gz"):
        base = base[:-3]
    if base.endswith(".fastq"):
        base = base[:-6]
    target = basedir + base + "_fastqc.html"
    return linkify(target, text=base)

### Parser for GTF files
### *** (this should be replaced with the parsers in genplot/genes)

def parseGTF(gtf):
    genes = {}
    transcripts = {}

    fn_open = gzip.open if gtf.endswith('.gz') else open
    with fn_open(gtf) as f:
        for line in f:
            if not line.startswith('#'):
                parseGTFrecord(line.rstrip("\r\n").split("\t"), genes, transcripts)
    return (genes, transcripts)

def parseGTFrecord(fields, genes, transcripts):
    what = fields[2]

    if what == 'gene':
        anndict = parseAnnotations(fields[8])
        gid = anndict['gene_id']
        genes[gid] = {'gene_id': gid,
                      'gene_biotype': anndict['gene_biotype'],
                      'gene_name': anndict['gene_name'] }

    elif what == 'transcript':
        anndict = parseAnnotations(fields[8])
        tid = anndict['transcript_id']
        transcripts[tid] = {'transcript_id': tid,
                            'gene_biotype': anndict['gene_biotype'],
                            'transcript_name': anndict['transcript_name'],
                            'gene_name': anndict['gene_name'] }

def parseAnnotations(ann):
    anndict = {}
    pieces = [ s.strip(" ") for s in ann.split(";") ]
    for p in pieces:
        pair = p.split(" ")
        if len(pair) == 2:
            anndict[pair[0]] = pair[1].strip('"')
    return anndict

### Adding annotations to a file

def annotateFile(infile, outfile, table, annot=[], annotNames=[], idcol=0, idname=False, missing='???'):
    """Read each line of file `infile' and look up the identifier in column `idcol' in 
dictionary `table'. If found, add the entries listed in `fields' to the current line.
Then write each line to `outfile'. If `idname' is specified, it is added to the front of the new header line
(to account for the fact that tables produced by R don't have a header for the first column)."""
    if len(annot) != len(annotNames):
        sys.stderr.write("annotateFile: the length of fields and fieldNames should be the same.\n")
        return
        
    m = [missing]*len(annot)
    print "Annotating file {} into {}.".format(infile, outfile)
    with open(outfile, "w") as out:
        with open(infile, "r") as f:
            hdr = parseLine(f.readline())
            if idname:
                hdr = [idname] + hdr
            ohdr = hdr[0:idcol+1] + annotNames + hdr[idcol+1:] # build output header line
            out.write("\t".join(ohdr) + "\n")
            for line in f:
                fields = parseLine(line)
                idv = fields[idcol].strip('"')
                if idv in table:
                    ann = table[idv]
                    ofields = fields[0:idcol+1]
                    for a in annot:
                        if a in ann:
                            ofields.append(ann[a])
                        else:
                            ofields.append(missing)
                    ofields += fields[idcol+1:]
                else:
                    ofields = fields[0:idcol+1] + m + fields[idcol+1:]
                out.write("\t".join(ofields) + "\n")

### For ma.py

def extractExpressions(sigfile, exprfile, outfile, sigidcol=0, sigfccol=1, idcol=0, datacol=1, maxrows=100):
    genes = []

    # Read list of genes with FCs
    with open(sigfile, "r") as f:
        f.readline()
        for line in f:
            parsed = line.rstrip("\r\n").split("\t")
            gid = parsed[sigidcol]
            fc = float(parsed[sigfccol])
            genes.append((abs(fc), gid))

    genes.sort(key=lambda g: g[0], reverse=True)
    wanted = [ g[1] for g in genes[:maxrows] ]

    with open(outfile, "w") as out:
        with open(exprfile, "r") as f:
            hdr = f.readline().rstrip("\r\n").split("\t")
            out.write("Gene\t" + "\t".join(hdr[datacol:]) + "\n")
            for line in f:
                parsed = line.rstrip("\r\n").split("\t")
                if parsed[0] in wanted:
                    out.write(parsed[idcol] + "\t" + "\t".join(parsed[datacol:]) + "\n")

### BAM utilities

def countReadsInBAM(filename):
    stats = pysam.idxstats(filename)
    nreads = 0
    for row in stats.split("\n"):
        row.rstrip("\r")
        fields = row.split("\t")
        if len(fields) > 2 and fields[0] != '*':
            nreads += int(fields[2])
    return nreads

def getChromsFromBAM(filename):
    chroms = []
    stats = pysam.idxstats(filename)
    for row in stats.split("\n"):
        fields = row.split("\t")
        if fields[0] != '*' and fields[0] != '':
            chroms.append(fields[0])
    return chroms

### Hub generator

class UCSCHub():
    genome = "hg38"
    name = "Hub"
    shortLabel = "testHub"
    longLabel = "A test UCSC hub"
    email = "ariva@ufl.edu"
    sizes = ""                  # Required by wigToBigWig, etc

    # These should not need to be changed
    dirname = "hub"             # Name of directory containing hub
    genomesFile = "genomes.txt"
    trackFile = "trackDb.txt"
    ucscpath = "/apps/dibig_tools/1.0/lib/ucsc/"
    parent = ""                 # If defined, we're in a container
    trackn = 1
    colorn = 0
    colors = []

    def __init__(self, name, shortLabel, longLabel, genome=None, sizes=None, dirname=None, genomesFile=None, trackFile=None, email=None):
        self.name = name
        self.shortLabel = shortLabel
        self.longLabel = longLabel
        self.trackn = 1
        self.initColors()
        if genome:
            self.genome = genome
        if dirname:
            self.dirname = dirname
        if trackFile:
            self.trackFile = trackFile
        if genomesFile:
            self.genomesFile = genomesFile
        if email:
            self.email = email
        if sizes:
            self.sizes = sizes
        else:
            self.sizes = self.ucscpath + "/" + self.genome + ".chrom.sizes"

    def initColors(self):
        self.colorn = 0
        self.colors = []
        for r in [255, 192, 64]:
            self.colors.append("{},0,0".format(r))
            self.colors.append("0,{},0".format(r))
            self.colors.append("0,0,{}".format(r))
            self.colors.append("{},{},0".format(r, r))
            self.colors.append("{},0,{}".format(r, r))
            self.colors.append("0,{},{}".format(r, r))

    def nextColor(self):
        c = self.colors[self.colorn]
        self.colorn += 1
        if self.colorn == len(self.colors):
            self.colorn = 0
        return c

    def mkdir(self, name):
        if not os.path.exists(name):
            os.mkdir(name)

    def missingOrStale(self, target, source):
        """Return True if `target' is missing or is older than `source'."""
        if os.path.isfile(target):
            if source and os.path.isfile(source):
                this = os.path.getmtime(target)
                that = os.path.getmtime(source)
                return (this < that)
            else:
                return False
        else:
            return True

    def shell(self, command, *args):
        """Like executef(), but allows multiple commands, redirection, piping, etc (runs a subshell).
Returns the command's output without the trailing \n."""
        cmd = command.format(*args)
        print "Executing: " + cmd
        try:
            return subprocess.check_output(cmd, shell=True).rstrip("\n")
        except subprocess.CalledProcessError as cpe:
            return cpe.output.rstrip("\n")

    def getBasename(self, path):
        return os.path.splitext(os.path.basename(path))[0]

    def generate(self):
        path = self.dirname + "/"                                  # hub/
        self.mkdir(path)
        self.mkdir(path + self.genome)                             # hub/hg38
        self.trackDbPath = self.genome + "/" + self.trackFile      # hg38/trackDb.txt

        # Write genomes file
        with open(path + self.genomesFile, "w") as out:            # hub/genomes.txt
            out.write("""genome {}
trackDb {}
""".format(self.genome, self.trackDbPath))

        # Write hub file
        with open(path + "hub.txt", "w") as out:                   # hub/hub.txt
            out.write("""hub {}
shortLabel {}
longLabel {}
genomesFile {}
email {}
""".format(self.name, self.shortLabel, self.longLabel, self.genomesFile, self.email))

        # Write trackDb file
        self.trackDbPath = path + self.trackDbPath
        with open(self.trackDbPath, "w") as out:                   # hub/hg38/trackDb.txt
            out.write("""# TrackDB for {} hub\n\n""".format(self.name))

    def writeStanza(self, out, base, options={}):
        if 'color' not in options:
            options['color'] = self.nextColor()
        out.write(base)
        for key, value in options.iteritems():
            out.write("{} {}\n".format(key, value))
        out.write("\n")

    def addWig(self, filename, name, shortLabel, longLabel, **options):
        bname = self.getBasename(filename)
        bw = bname + ".bw"
        bwpath = self.dirname + "/" + self.genome + "/" + bw
        if self.missingOrStale(bwpath, filename):
            self.shell("{}/wigToBigWig -clip {} {} {}", self.ucscpath, filename, self.sizes, bwpath)

        with open(self.trackDbPath, "a") as out:
            self.writeStanza(out, """track {}
{}bigDataUrl {}
shortLabel {}
longLabel {}
type bigWig
autoScale on
priority {}
""".format(name, self.parent, bw, shortLabel.format(name), longLabel.format(name), self.trackn), options=options)
            self.trackn += 1

    def addBedGraph(self, filename, name, shortLabel, longLabel, **options):
        bname = self.getBasename(filename)
        bw = name + ".bw"
        bwpath = self.dirname + "/" + self.genome + "/" + bw
        if self.missingOrStale(bwpath, filename):
            print self.shell("{}/bedGraphToBigWig {} {} {}", self.ucscpath, filename, self.sizes, bwpath)

        with open(self.trackDbPath, "a") as out:
            self.writeStanza(out, """track {}
{}bigDataUrl {}
shortLabel {}
longLabel {}
type bigWig
autoScale on
priority {}
""".format(name, self.parent, bw, shortLabel.format(name), longLabel.format(name), self.trackn), options=options)
            self.trackn += 1

    def addBed(self, filename, name, shortLabel, longLabel, **options):
        bname = self.getBasename(filename)
        bw = bname + ".bw"
        bwpath = self.dirname + "/" + self.genome + "/" + bw
        if self.missingOrStale(bwpath, filename):
            self.shell("{}/bedToBigBed -clip {} {} {}", self.ucscpath, filename, self.sizes, bwpath)

        with open(self.trackDbPath, "a") as out:
            self.writeStanza(out, """track {}
{}bigDataUrl {}
shortLabel {}
longLabel {}
type bigWig
autoScale on
priority {}
""".format(name, self.parent, bw, shortLabel, longLabel, self.trackn), options=options)
            self.trackn += 1

    def startContainer(self, name, shortLabel, longLabel, tracktype):
        with open(self.trackDbPath, "a") as out:
            out.write("""track {}
compositeTrack on
longLabel {}
shortLabel {}
type {}
allButtonPair on

""".format(name, longLabel, shortLabel, tracktype))
        self.parent = "parent {} on\n".format(name)

    def endContainer(self):
        self.parent = ""
