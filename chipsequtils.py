# (c) 2017, A. Riva, DiBiG
# RNA-Seq Analysis pipeline helper functions, v2.0

import Utils

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
