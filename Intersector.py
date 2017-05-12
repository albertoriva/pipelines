### This class generates all pairwise intersections
### between up- and down-regulated entries in a set
### of contrasts

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
                interup = "{}and{}-up.csv".format(a, b)
                interdn = "{}and{}-dn.csv".format(a, b)
                out.write("colx.py -o {} {} {}\n".format(interup, da['up'], db['up']))
                out.write("colx.py -o {} {} {}\n".format(interdn, da['down'], db['down']))

