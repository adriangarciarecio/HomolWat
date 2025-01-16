import MySQLdb as mdb
import glob, os



# entra en pdblib para coger la info

try:
    con = mdb.connect('alf03.uab.cat', 'lmcdbuser', 'lmcdbuser12345', 'lmcdb')
except:
    print 'could not connect'
    exit(0)

with con:
    cur = con.cursor()
    cur.execute("SELECT pdbid, uniprot, domain, inactive_set, filt_wat, total_wat, resolution, date, chain, taxo, state FROM pdb")
    rows = cur.fetchall()

# guarda la info que nos interesa en un diccionario (Rec_pdb_dict)
recs, pdbs = [], []
Rec_pdb_dict = {}
cont = 0


Class_GPCR = {'7tm_1':'A', '7tm_2':'B', '7tm_3':'C', 'Frizzled':'F'}
##print "num pdbid receptor inn_wats wats resolution"
for x in rows:
    pdbid = x[0]
    uniprot = x[1]
    domain = x[2]
    inner_wats = x[4]
    l_wats = inner_wats.split(",")
    inact_set = x[3]
    water = x[5]
    l_all_wats = water.split(",")
    resolution = x[6]
    date = x[7]
    chain = x[8]
    chains = chain.replace(" ", "")
    l_chains = chains.split(",")
    org = x[9]
    state = x[10]
    recep = x[1].split('_')
    if domain in ['7tm_1', '7tm_2', '7tm_3', 'Frizzled']:
        for i, numW in enumerate(l_wats):
            if int(numW) > 0:
                print pdbid, recep[0], recep[1], Class_GPCR[domain], numW, l_all_wats[i], resolution, date, state, org
    # if domain == '7tm_1':
        #if inact_set == 1:
        # if 1 == 1:
    
        #if recep[0] not in recs:
         #   recs.append(recep[0])
        #    pdbs.append(x[0])
        # cont += 1
#        if resolution <= 2.5:
	# if inner_wats > 0:
	# print pdbid, recep[0], recep[1], inner_wats, water, resolution, date, state, org

