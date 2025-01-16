###############################################################################
'''

This script UPDATE the database of crystals in /REC_WAT/ and files for the
blast of model vs GPCRdb to assign subfamily to model
files: ref_gpcr.phr, ref_gpcr.pin and ref_gpcr.psq

Info of crystals is taken from LMCDB
Eduardo Mayol (e.mayoles@gmail.com)
Jul19
'''
###############################################################################


path_lib = "/people/common/add_water/REC_WAT/"
path_in = "/people/common/LMCDB/PDB/REC_WAT/"
path_aln = "/people/common/add_water/aln/"
# path_alf = "/var/www/lmc.uab.es/HW/hw/tmp/REC_WAT/"

import MySQLdb as mdb
import shutil

try:
    con = mdb.connect('alf03.uab.cat', 'lmcdbuser', 'lmcdbuser12345', 'lmcdb')
except:
    print('could not connect')
    exit(0)

with con:
    cur = con.cursor()
    cur.execute("SELECT pdbid, chain, domain, inactive_set, filt_wat, uniprot FROM pdb")
    rows = cur.fetchall()

pdb_list = []
codes = []
cont, all_wats = 0, 0
for row in rows:
    pdbid = row[0]
    chain = row[1]
    chains = chain.replace(" ", "")
    l_chains = chains.split(",")
    dom = row[2]
    inact = row[3]
    filt_wats = row[4]
    l_wats = filt_wats.split(",")
    uniprot = row[5]
    unip = uniprot.split("_")
    family = unip[0]
    # print(pdbid, l_wats, len(l_wats), l_chains, len(l_chains))
    if dom in ['7tm_1', '7tm_2', '7tm_3', 'Frizzled']: 
        for i, numW in enumerate(l_wats):
            if int(numW) > 0:
                code = pdbid + "_" + l_chains[i] + "_" + family 
                codes.append(code)
                print(code , numW, dom)
                cont += 1
                all_wats += int(numW)
print("ch", cont, "wats", all_wats)
#copy pdbs with waters
for pdb in codes:
    shutil.copyfile(path_in + pdb + ".pdb", path_lib + pdb + ".pdb")


