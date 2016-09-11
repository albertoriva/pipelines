# (c) 2015, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import sys
import Library

REGISTRY = {'test':       Library.testLine,
            'samples':    Library.Samples,
            'rnasamples': Library.RNAseqSamples,
            'trim':       Library.Trimmer,
            'fastqcount': Library.FASTQcounter,
            'star':       Library.StarAligner,
            'startx':     Library.StarAlignerTx,
            'bamcount':   Library.BAMcounter,
            'markdup':    Library.Markdup,
            'merge':      Library.BAMmerger,
            'bamcat':     Library.BAMconcatenator,
            'rsemquant':  Library.RSEMquant,
            'rsemdiff':   Library.RSEMdiff,
            'bamtowig':   Library.BAMtoWig,
            'macs':       Library.MACScallpeak,
            'csfilter':   Library.CSfilter,
            'mmap':       Library.MMAP,
            'cscall':     Library.CScall,
            'mcomp':      Library.MCOMP,
            'genemeth':   Library.GeneMeth,
            'methplots':  Library.MethPlots,
            'rencontig':  Library.ContigRenamer,
            'refindex':   Library.IndexReference,
            'prokka':     Library.Prokka,
            'mauve':      Library.Mauve,
            'progmauve':  Library.ProgMauve,
            'roary':      Library.Roary,
            'snpsites':   Library.SNPsites
}

class Director():
    """This class coordinates the execution of the pipeline."""

    actor = None
    steps = []

    def __init__(self, actor):
        self.actor = actor

    def add(self, key, **properties):
        dkey = key
        p = key.find(".")
        if p > -1:
            dkey = key[:p]

        if dkey in REGISTRY:
            cls = REGISTRY[dkey]
            line = cls(self.actor, key=key, properties=properties)
            self.steps.append(line)
        else:
            sys,stderr.write("Warning: no Line with key `{}'.\n".format(dkey))
        return line

    def startAt(self, startkey):
        """Set all steps in this pipeline to dry, until `startkey' is reached. All steps
from `startkey' onwards will be set to not-dry."""
        if startkey:
            dry = True
            for s in self.steps:
                if s.key == startkey:
                    dry = False
                s.dry = dry

    def stopAt(self, stopkey):
        """Set all steps after `stopkey' to dry."""
        print "Stopping at {}".format(stopkey)
        if stopkey:
            ok = False
            for s in self.steps:
                # print s.key
                if ok:
                    # print "Setting {} to dry".format(s.name)
                    s.dry = True
                elif s.key == stopkey:
                    # print "Found {}".format(s.key)
                    ok = True

    def showSteps(self, wait=True):
        print "Ready to run the following steps:"
        for s in self.steps:
            print "{} {}".format("-" if s.dry else "+", s.name)
        if wait:
            print "Press Enter to start execution."
            raw_input()

    def PerformAll(self, method, immediatestop=False):
        good = True
        for l in self.steps:
            self.actor.log.log("Director: performing {} on `{}'.", method, l.name)
            m = getattr(l, method)
            f = m()
            if not f:
                sys.stderr.write("Error in {}: {}: {}\n".format(method, l.name, l.status))
                if immediatestop:
                    return False
                else:
                    good = False
        return good

    def VerifyAll(self):
        return self.PerformAll('Verify')

    def PreExecuteAll(self):
        return self.PerformAll('PreExecute')

    def ExecuteAll(self):
        return self.PerformAll('Execute', immediatestop=True)

    def PostExecuteAll(self):
        return self.PerformAll('PostExecute')

    def ReportAll(self):
        return self.PerformAll('Report')

    def RunScript(self):
        if not self.VerifyAll():
            return False
        if not self.PreExecuteAll():
            return False
        if not self.ExecuteAll():
            return False
        if not self.PostExecuteAll():
            return False
        if not self.ReportAll():
            return False

