# (c) 2015, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import os

class Line():
    """This class is the parent of all classes representing operations performed by the pipelines.
Its methods all do nothing and return True."""
    
    name = "(no name)"
    key = ""
    actor = None
    dry = False
    status = ""
    waiters = []
    properties = {}
    tempfiles = []

    def __init__(self, act, key="", properties={}):
        self.actor = act
        self.key = key
        self.status = ""
        self.waiters = []
        self.tempfiles = []
        self.properties = properties
        if 'dry' in properties:
            self.dry = properties['dry']
        self.Setup()

    def tempfile(self, filename):
        self.tempfiles.append(filename)
        return filename

    def cleanTempfiles(self):
        for f in self.tempfiles:
            try:
                os.remove(f)
            except:
                pass

    def Setup(self):
        """The Setup() method is called by init."""
        return True

    def Verify(self):
        """The Verify() method is called before starting the pipeline. It should check
for prerequisites (e.g., necessary software tools) and return False if they are not 
met. This prevents the pipeline from getting started if we already know it'll fail
down the road."""
        return True

    def PreExecute(self):
        """The PreExecute() method is called when the pipeline is ready to start. The PreExecute
methods of all Lines in the Script are called in order. This method can be used for setup 
operations (e.g., letting subsequent Lines know that this one will be executed)."""
        return True

    def Execute(self):
        """The Execute() method performs the actions of this Line. It should honor the 'rehearse'
flag (if true, do everything except actually running the actions)."""
        return True

    def PostExecute(self):
        """The PostExecute method should take care of cleanup, checking for job completion, etc. 
All PostExecute methods are called in order after the pipeline has finished executing."""
        return True

    def Report(self):
        """The Report method writes the report to `stream'."""
        return True


        
