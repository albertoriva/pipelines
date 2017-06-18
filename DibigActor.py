###################################################
#
# (c) 2016, Alberto Riva, ariva@ufl.edu
# DiBiG, ICBR Bioinformatics, University of Florida
#
# See the LICENSE file for license information.
###################################################

from MultiSampleActor import MultiSampleActor
from Logger import Logger

class DibigActor(MultiSampleActor):

    def __init__(self):
        self._addFile(self.libpath + "js/jquery.tablescroll.js") # *** The first two should be inherited!
        self._addFile(self.libpath + "css/jquery.tablescroll.css")
        self._addFile(self.libpath + "img/UF-ICBR-logo.png")

    def initFiles(self):
        self.log.setLogfile(self.getConf("logfile"))
        self.log.setEcho('stdout')
        self.log.logStart(self.title)

        ## Ensure we don't have old .done files and tmp- files lying around
        self.shell("rm -f *.done tmp-* .files")

        ## Initialize .files
        self._addToInclude("*.html", "*.png", "*.pdf", "*.xlsx", "*.csv", "*.css", "*.js", "*.bed", "*.vcf", "*.bedGraph", "*.conf")

    def headExtra(self):
        """Returns additional tags for the <HEAD> section."""
        return """    <link rel="stylesheet" type="text/css" href="jquery.tablescroll.css"/>
    <!-- <link rel="stylesheet" type="text/css" href="bootstrap.min.css"/> -->
    <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.5.1/jquery.min.js"></script>
    <script type="text/javascript" src="jquery.tablescroll.js"></script>
"""  

    def header(self):
        """Returns the header part of the HTML report (called by preamble)."""
        return """<table class='hdr'>
      <tr><td align='left'><A href='http://biotech.ufl.edu/'><img src='UF-ICBR-logo.png' border='0'></A></td><td align='right'><A class='dibig' href='http://dibig.biotech.ufl.edu'>DiBiG</A></td></tr>
      <tr><th class='hdr' align='left'>ICBR Bioinformatics</th><th class='hdr' align='right'><i>Powered by Actor, v1.0</i><tr>
    </table>"""

