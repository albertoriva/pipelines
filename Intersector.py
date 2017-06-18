### This class generates all pairwise intersections
### between up- and down-regulated entries in a set
### of contrasts

import Table

def triangleGenerator(keys):
    nk = len(keys)
    for i1 in range(nk):
        for i2 in range(i1+1, nk):
            yield (keys[i1], keys[i2])

class Intersector():
    contrasts = []
    keys = []
    data = {}
    label = ""
    geneCol = 2
    fcCol = 3

    outfiles = {}

    def __init__(self, contrasts, label, geneCol=2, fcCol=3):
        """Initialize this intersector for the specified `contrasts'. `label'
is one of genefinal, codingfinal, isofinal."""
        self.contrasts = contrasts
        self.keys = []
        self.data = {}
        self.triangle = {}
        self.label = label
        self.geneCol = geneCol
        self.fcCol = fcCol
        self.outfiles = {}

        idx = 65
        for contr in contrasts:
            key = chr(idx)
            self.keys.append(key)
            vs = "{}.vs.{}".format(contr['test'], contr['control'])
            self.data[key] = {'vs': vs,
                              'source': contr[self.label],
                              'up': vs + "-up.csv",
                              'down': vs + "-down.csv"}
            idx += 1

        # Initialize triangle
        for ai in range(len(self.keys)):
            adict = {}
            self.triangle[self.keys[ai]] = adict
            for bi in range(ai+1, len(self.keys)):
                adict[self.keys[bi]] = {'key': "{}and{}".format(self.keys[ai], self.keys[bi])}

        # print self.keys
        # print self.data
        # print self.triangle
        # raw_input()

    def splitUpDown(self):
        for k in self.keys:
            d = self.data[k]
            nup = 0
            ndown = 0
            with open(d['up'], "w") as fup:
                with open(d['down'], "w") as fdown:
                    with open(d['source'], "r") as f:
                        f.readline() # Skip header
                        for line in f:
                            parsed = line.rstrip("\r\n").split("\t")
                            if float(parsed[self.fcCol]) > 0:
                                fup.write(parsed[self.geneCol] + "\n")
                                nup += 1
                            else:
                                fdown.write(parsed[self.geneCol] + "\n")
                                ndown += 1
            d['nup'] = nup
            d['ndown'] = ndown

    def writeIntersectScript(self, filename):
        with open(filename, "w") as out:
            out.write("""#!/bin/bash\n\nmodule load dibig_tools\n\n""")
            for (a, b) in triangleGenerator(self.keys):
                da = self.data[a]
                db = self.data[b]
                upA = "{}{}-{}only-up.csv".format(a, b, a)
                upB = "{}{}-{}only-up.csv".format(a, b, b)
                up0 = "{}{}-{}and{}-up.csv".format(a, b, a, b)
                dnA = "{}{}-{}only-dn.csv".format(a, b, a)
                dnB = "{}{}-{}only-dn.csv".format(a, b, b)
                dn0 = "{}{}-{}and{}-dn.csv".format(a, b, a, b)
                out.write("colx.py -o {} {} {}\n".format(up0, da['up'], db['up']))
                out.write("colx.py -o {} -d {} {}\n".format(upA, da['up'], db['up']))
                out.write("colx.py -o {} -d {} {}\n".format(upB, db['up'], da['up']))
                out.write("colx.py -o {} {} {}\n".format(dn0, da['down'], db['down']))
                out.write("colx.py -o {} -d {} {}\n".format(dnA, da['down'], db['down']))
                out.write("colx.py -o {} -d {} {}\n".format(dnB, db['down'], da['down']))
                self.outfiles[(a, b)] = (upA, up0, upB, dnA, dn0, dnB)

    def toHTML(self, ACT, tblid):
        tblid1 = tblid + "_1"
        tbl1 = Table.ScrollingTable(id=tblid1, align="LLRC", caption="Intersections of over-expressed genes.")
        tbl1.startHead()
        tbl1.addHeaderRow(["Contrast 1", "Contrast 2", "1 only - N", "1 only - genes", "Common - N", "Common - genes", "2 only - N", "2 only - genes"])
        tbl1.startBody()
        for pair in triangleGenerator(self.keys):
            (a, b) = pair
            files = self.outfiles[pair]
            nupA = ACT.fileLines(files[0])
            nup0 = ACT.fileLines(files[1])
            nupB = ACT.fileLines(files[2])

            tbl1.addRow([ Table.C(a + " - " + self.data[a]['vs']), Table.C(b + " - " + self.data[b]['vs']),
                          Table.D(nupA, cls='upreg'),
                          Table.A(files[0]),
                          Table.D(nup0, cls='upreg'),
                          Table.A(files[1]),
                          Table.D(nupB, cls='upreg'),
                          Table.A(files[2]) ])

        tbl1.toHTML(ACT.out)

        tblid2 = tblid + "_2"
        tbl2 = Table.ScrollingTable(id=tblid2, align="LLRC", caption="Intersections of under-expressed genes.")
        tbl2.startHead()
        tbl2.addHeaderRow(["Contrast 1", "Contrast 2", "1 only - N", "1 only - genes", "Common - N", "Common - genes", "2 only - N", "2 only - genes"])
        tbl2.startBody()
        for pair in triangleGenerator(self.keys):
            (a, b) = pair
            files = self.outfiles[pair]
            nupA = ACT.fileLines(files[3])
            nup0 = ACT.fileLines(files[4])
            nupB = ACT.fileLines(files[5])

            tbl2.addRow([ Table.C(a + " - " + self.data[a]['vs']), Table.C(b + " - " + self.data[b]['vs']), 
                          Table.D(nupA, cls='dnreg'),
                          Table.A(files[3]),
                          Table.D(nup0, cls='dnreg'),
                          Table.A(files[4]),
                          Table.D(nupB, cls='dnreg'),
                          Table.A(files[5]) ])

        tbl2.toHTML(ACT.out)
