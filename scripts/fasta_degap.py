"""
open a fasta files with gaps(-) and delete them
"""

import MySQLdb as mdb

# entra en pdblib para coger la info

try:
    con = mdb.connect('alf03.uab.cat', 'lmcdbuser', 'lmcdbuser12345', 'lmcdb')
except:
    print 'could not connect'
    exit(0)

with con:
    cur = con.cursor()
    cur.execute("SELECT pdbid, domain, inactive_set FROM pdb")
    rows = cur.fetchall()

pdbs = []
for x in rows:
    if x[1] == '7tm_1':
        if x[2] == 1:
            pdbs.append(x[0])

## print list of ref pdbs
print (pdbs, len(pdbs))


path = "/home/edu/Documentos/LCM/water_molec/align/"
fasta_wgap = "xrays_2fasta.txt"

f_in = open(path + fasta_wgap, "r")
f_r = f_in.readlines()
f_in.close()

fasta_out = []
rec = ""
cont = 0
c1 = 0
line_out = ""
for line in f_r:
    if line.startswith(">"):
        c1 += 1
        full_code = line.strip()
        code = full_code.split("_")
        pdb = code[-1]
        receptor = code[0][1:]
        print line[:-1]
        #print pdb, receptor
        if line_out != "":
            fasta_out.append(line_out)
            line_out = ""
        if pdb not in pdbs:
            info = ">" + receptor + "_DELETE"
            fasta_out.append(info)
        else:
            info = ">" + receptor
            fasta_out.append(info)

    else:
        old_line = line.strip()
        #print old_line
        for let in old_line:
            #print let
            if let != "-":
                line_out += let
                #print line_out, len(line_out)
                if len(line_out) == 60:
                    fasta_out.append(line_out)
                    line_out = ""
print c1, "recs in fasta"
print "NEEEW!!!!"
f_out = open(path + "ref_xrays.fasta", "w")
for l_o in fasta_out:
    if l_o != fasta_out[-1]:
        f_out.write(l_o + "\n")
    else:
        f_out.write(l_o)
f_out.close()