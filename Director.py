# (c) 2015, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import sys
import importlib

class Director():
    """This class coordinates the execution of the pipeline."""

    actor = None
    steplist = []               # Names of steps
    steps = []                  # Actual step objects
    registry = {}

    def __init__(self, actor, library="Library"):
        self.actor = actor
        self.steplist = []
        self.steps = []

        lib = importlib.import_module(library)
        self.registry = lib.REGISTRY

    def setSteps(self, steplist):
        """Set the list of steps to be performed by this dirctor to `steplist'. Steplist
can be either a list of strings or a string containing comma-separated step names (e.g.
"step1, step2, step3". Use the step() method to know if a step should be executed."""
        if type(steplist).__name__ == 'str':
            steplist = [ i.strip(" ") for i in steplist.split(",") ]
        self.steplist = steplist
        self.notifiedSteps = []

    def stepPresent(self, step):
        return (step in self.steplist) or ("-"+step in self.steplist) or ("no"+step in self.steplist)

    def stepDry(self, step):
        return ("-"+step in self.steplist) or ("no"+step in self.steplist)

    # def step(self, wanted):
    #     if wanted in self.steplist:
    #         if not wanted in self.notifiedSteps: # should we notify?
    #             print "Performing step `{}'.".format(wanted)
    #             self.notifiedSteps.append(wanted)
    #         return True
    #     else:
    #         if not wanted in self.notifiedSteps:
    #             print "Skipping step `{}'.".format(wanted)
    #             self.notifiedSteps.append(wanted)
    #         return False

    def step(self, key, **properties):
        if self.stepPresent(key):
            self.add(key, dry=self.stepDry(key), **properties)

    def add(self, key, **properties):
        line = None
        dkey = key
        p = key.find(".")
        if p > -1:
            dkey = key[:p]

        if dkey in self.registry:
            cls = self.registry[dkey]
            line = cls(self.actor, key=key, properties=properties)
            self.steps.append(line)
        else:
            sys.stderr.write("Warning: no Line with key `{}'.\n".format(dkey))
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

    def showSteps(self):
        print "Ready to run the following steps:"
        for s in self.steps:
            print "{} {}".format("-" if s.dry else "+", s.name)
        if self.actor.ask:
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

