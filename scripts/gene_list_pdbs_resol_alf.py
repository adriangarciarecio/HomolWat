import MySQLdb as mdb
import glob,os

try:
    con = mdb.connect('alf03.uab.cat', 'lmcdbuser', 'lmcdbuser12345', 'lmcdb')
except:
    print 'could not connect'
    exit(0)

with con:
    cur = con.cursor()
    cur.execute("SELECT pdbid, uniprot, domain, inactive_set, chain, resolution, filt_wat, state FROM pdb")
    rows = cur.fetchall()

# guarda la info que nos interesa en un diccionario (Rec_pdb_dict)
recs, pdbs, chains = [], [], []
Rec_pdb_dict_ACT, Rec_pdb_dict_INACT = {}, {}
for x in rows:
    if x[2] in ['7tm_1', '7tm_2', '7tm_3', 'Frizzled']:
        if x[3] == 1:
            recep = x[1].split('_')
            if recep[0] not in recs:
                recs.append(recep[0])
                pdbs.append(x[0])

#print sorted(recs)
#print sorted(pdbs)

c = 0
for r in sorted(recs):
    #print '---', r, '---'
    sort_resol_act,sort_resol_inact, pdbs_act, pdbs_inact = [], [], [], []
    for y in rows:
        rec = y[1].split('_')
        if rec[0] == r:
            c += 1
            inner_wats = y[6]
            l_wats = inner_wats.split(",")
            for i, numW in enumerate(l_wats):
                if int(numW) > 0:  # internal wats
                    info = [y[0], y[5]]  ## pdb, resol, state
                    if y[7].startswith('inact'):
                        if info not in sort_resol_inact:
                            sort_resol_inact.append(info)
                    else:
                        if info not in sort_resol_act:
                            sort_resol_act.append(info)
            #pdbs.append(y[0])
    sorted_l_act = sorted(sort_resol_act, key=lambda resol: resol[1])
    sorted_l_inact = sorted(sort_resol_inact, key=lambda resol: resol[1])
    for elem in sorted_l_act:
        pdbs_act.append(elem[0])
    for elem in sorted_l_inact:
        pdbs_inact.append(elem[0])
    Rec_pdb_dict_ACT[r] = pdbs_act + pdbs_inact
    Rec_pdb_dict_INACT[r] = pdbs_inact + pdbs_act
    #print r
#print len(recs)

#path = "/people/common/add_water/"
###print file_exist('2Y00_A')
##folder = glob.glob(path + "REC_WAT/" + "*.pdb")

#sodiums = []
f_outAct = open("Pdbs_resol_ACT", "w")
f_outInact = open("Pdbs_resol_INACT", "w")
recs_resol_act, recs_resol_inact = [], []

for k, v in Rec_pdb_dict_ACT.items():
    ## k is for family code, v for pdb_code LIST
    info = str(k + "_" + ' '.join(v))
    print info
    recs_resol_act.append(info)

for k, v in Rec_pdb_dict_INACT.items():
    ## k is for family code, v for pdb_code LIST
    info = str(k + "_" + ' '.join(v))
    print info
    recs_resol_inact.append(info)

for line in recs_resol_act:
    if line != recs_resol_act[-1]:
        f_outAct.write(line + "\n")
    else:
        f_outAct.write(line)
f_outAct.close()

for line in recs_resol_inact:
    if line != recs_resol_inact[-1]:
        f_outInact.write(line + "\n")
    else:
        f_outInact.write(line)
f_outInact.close()