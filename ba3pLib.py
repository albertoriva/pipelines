# (c) 2015, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import os

import Utils
from SampleCollection import SampleCollection
from Lines import Line
import Table

### BA3P

### Samples manager

class Samples(Line):
    """This object handles initializing the SampleCollection, checking that all files exist, fixing paths."""
    name = "SamplesManager"

    def Setup(self):
        ACT = self.actor
        ACT.sc = SampleCollection(ACT.Conf)
        SC = ACT.sc

        ## Fix paths to fastq files in each readset and store base name.
        #for rs in SC.readsets:
        #    rs['left'] = ACT.checkPath(rs['left'])
        #    rs['lbase'] = ACT.setFileExt(os.path.split(rs['left'])[1], "", remove=[".fastq", ".gz"]) # basename of left-side reads
        #    if rs['paired']:
        #        rs['right'] = ACT.checkPath(rs['right'])
        #        rs['rbase'] = ACT.setFileExt(os.path.split(rs['right'])[1], "", remove=[".fastq", ".gz"]) # basename of right-side reads

        ## Get fasta names from conf file for samples
        for smp in SC.samples:
            smp['fasta'] = "../" + ACT.getConf("fasta", section=smp['name'])

        return True

    def Verify(self):
        return self.actor.sc.verify()

    def Report(self):
        ACT = self.actor
        SC = ACT.sc

        ACT.scene("Input data")
        ACT.reportf("The following table summarizes the samples, conditions, and contrasts in this analysis. A readset is either a single fastq file or a pair of fastq files (for paired-end sequencing).")
        tbl1 = Table.ScrollingTable(id='tbl1', caption='Summary of input data')
        tbl1.startHead()
        tbl1.addHeaderRow(["Category", "Data"])
        tbl1.startBody()
        tbl1.addSectionRow("Summary of input data")
        tbl1.addRowHeader("Contrasts:")
        tbl1.addRow([ ", ".join([ contr['test'] + " vs. " + contr['control'] for contr in SC.contrasts])])
        tbl1.addRowHeader("Conditions:")
        tbl1.addRow([ ", ".join([ cond['name'] for cond in SC.conditions])])
        for cond in SC.conditions:
            condsamples = SC.conditionSamples(cond)
            tbl1.addRowHeader(cond['name'], rowspan=len(condsamples))
            for smp in condsamples:
                tbl1.addRow([ "{} ({} readsets)".format(smp['name'], len(smp['readsets'])) ])
        tbl1.toHTML(ACT.out)
        return True

class ContigRenamer(Line):
    """Rename sequences in contig files as XXX_n, where XXX is the sample name and n is a progressive number."""
    name = "Rename contigs"

    def Execute(self):
        ACT = self.actor
        SC = ACT.sc

        ACT.log.log("Creating directory Contigs/")
        ACT.mkdir("Contigs")

        ACT.shell("touch dummymap")
        for smp in SC.samples:
            infasta = smp['fasta']
            outfasta = "Contigs/" + smp['name'] + ".spades.fasta"
            print outfasta
            smp['contig'] = outfasta
            if not self.dry:
                ACT.shell("module load dibig_tools; fastools -map -mapfile-in dummymap -in {} -out {} -if-missing counter -prefix {}_", infasta, outfasta, smp['name'])
        return True

class IndexReference(Line):
    """This line checks if the reference file in ACT.reference needs to be indexed with samtools."""
    name = "Index reference"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log

        fai = ACT.reference + ".fai" 
        refdict = ACT.setFileExt(ACT.reference, ".dict")

        if ACT.missingOrStale(fai, ACT.reference):
            LOG.log("Reference index {} does not exist or is out of date, creating it.", fai)
            if not self.dry:
                ACT.shell("module load samtools; samtools faidx {}", ACT.reference)

        if ACT.missingOrStale(refdict, ACT.reference):
            LOG.log("Reference dictionary {} does not exist, creating it.", refdict)
            if not self.dry:
                ACT.submit("picard.qsub CreateSequenceDictionary R={} O={}".format(ACT.reference, refdict), done="dict.done")
                ACT.wait("dict.done")
        return True

class Prokka(Line):
    name = "Prokka"

    def Verify(self):
        ACT = self.actor

        if ACT.genus == None or ACT.species == None:
            self.status = "prokka requires specification of genus and species in configuration file."
            return False
        else:
            return True
        
    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        ACT.nprokka = 0
        for smp in SC.samples:
            prokkadir = smp['name']
            smp['prokkadir'] = prokkadir
            fasta = smp['name'] + ".spades.fasta"
            if ACT.locus == None:
                locus = ""
            else:
                locus = "locus=" + ACT.locus
            if not self.dry:
                os.chdir("Contigs")
                ACT.submit('prokka.qsub {} {} "{}" "{}" {} prefix={}'.format(fasta, prokkadir, ACT.genus, ACT.species, locus, smp['name']), done="../prokka.@.done")
                os.chdir("..")
                ACT.nprokka += 1
        return True

class Mauve(Line):
    """Run mauve contig mover on all fasta files in the Contigs directory."""
    name = "Mauve contig mover"

    def Execute(self):
        ACT = self.actor
        LOG = ACT.log
        SC = ACT.sc

        LOG.log("Running Mauve contig mover on directory Contigs/")
        nmauve = 0
        if not self.dry:
            # ACT.submit("mauve-contig-mover.qsub {}".format(ACT.reference), done="mauve.done")
            # ACT.wait("mauve.done")
            for smp in SC.samples:
                ACT.submit("mauve-contig-mover.qsub {} {}".format(ACT.reference, smp['contig']), done="mauve.@.done")
                nmauve += 1
            ACT.wait(("mauve.@.done", nmauve))
        LOG.log("Selecting optimal alignment for each sample.")
        ACT.mkdir("Alignments")
        for smp in SC.samples:
            base = "Ordered/" + smp['name'] + ".spades/"
            LOG.log("{}: examining directory {}.", smp['name'], base)
            alns = os.listdir(base)
            maxa = [0, '']
            for a in alns:
                if a[0:9] == 'alignment':
                    i = int(a[9:])
                    if i > maxa[0]:
                        maxa = [i, a]
            LOG.log("{}: best alignment is {}", smp['name'], maxa[0])
            ACT.shell("cp {}/{}/{}.spades.fasta Alignments/".format(base, maxa[1], smp['name']))

        return True
    
class ProgMauve(Line):
    """Run progressive-mauve on all files in the Alignments directory. Files listed in the 'additional'
property are added to the directory before running progressive mauve."""
    name = "Progressive Mauve"
    
    def Execute(self):
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log
        
        if not self.dry:
            os.chdir("Alignments")
            if ACT.additional:
                for a in ACT.additional:
                    if os.path.exists(a):
                        LOG.log("Copying additional reference {} into Alignments.", a)
                        ACT.shell("cp {} .".format(a))
                    else:
                        LOG.log("Additional reference {} not found.", a)
            if ACT.mauvemem:
                ACT.submit("progressive-mauve.qsub denovo-final denovo-final", done="../progr.done", options="--time="+ACT.mauvetime)
            else:
                ACT.submit("progressive-mauve.qsub denovo-final denovo-final", done="../progr.done")
            os.chdir("..")
            ACT.wait("progr.done")
        return True

    def PostExecute(self):
        ACT = self.actor
        
        ## Convert final xmfa to fasta
        ACT.shell("module load dibig_tools perl; convertAlignment.pl -i Alignments/denovo-final.xmfa -o denovo-final.fasta -f fasta -g xmfa -c")
        return True

class Roary(Line):
    name = "Roary"
    
    def Execute(self):
        ACT = self.actor
        SC = ACT.sc
        LOG = ACT.log

        ACT.wait(("prokka.@.done", ACT.nprokka))
        if True: # ACT.nprokka > 0:
            ACT.mkdir("Roary")
            for smp in SC.samples:
                ACT.shell("cp Contigs/{}/*.gff Roary/".format(smp['prokkadir']))
            if not self.dry:
                os.chdir("Roary")
                ACT.submit("roary.qsub", done="../roary.done")
                os.chdir("..")
                ACT.wait("roary.done")

            ACT.shell("cp Roary/core_gene_alignment.aln .")
            ACT.shell("module load snp-sites; snp-sites -m -o core-alignment-snps.fa core_gene_alignment.aln")
        return True

class SNPsites(Line):
    name = "snp-sites"
    infile = None
    outfile = None

    def Setup(self):
        self.infile = Utils.dget('infile', self.properties, self.infile)
        self.outfile = Utils.dget('outfile', self.properties, self.outfile)

    def Verify(self):
        if self.infile and self.outfile:
            return True
        self.status = "Need to specify infile and outfile!"
        return False

    def Execute(self):
        ACT = self.actor

        ACT.shell("module load snp-sites; snp-sites -m -o {} {}".format(self.outfile, self.infile))
        return True

REGISTRY = {'samples':    Samples,
            'rencontig':  ContigRenamer,
            'refindex':   IndexReference,
            'prokka':     Prokka,
            'mauve':      Mauve,
            'progmauve':  ProgMauve,
            'roary':      Roary,
            'snpsites':   SNPsites
}
