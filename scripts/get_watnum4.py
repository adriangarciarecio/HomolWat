#!/usr/bin/env python

"""
check the number of waters and internal waters of pdbs
"""

import glob
import numpy as np
import sys
import os
import pandas as pd

path = os.path.dirname(os.path.abspath("__file__"))
sys.path.insert(0, path)

def check_CircVar(CB_coords, wat_list):
    cont = 0
    rcv = 10.0
    for w in wat_list:
        coxB = w[30:38].strip()
        coyB = w[38:46].strip()
        cozB = w[47:54].strip()
        coord_w = np.array([coxB, coyB, cozB], float)
        sumvec = 0
        sumcount = 0
        for cb in CB_coords:
            coord_cb = np.array(cb, float)
            vec = coord_w - coord_cb
            modvec = np.linalg.norm(vec)
            if modvec < rcv:
                sumvec += vec / modvec
                sumcount += 1
        cv = 1 - (np.linalg.norm(sumvec)) / sumcount
        #print w_num, cv
        if cv > 0.6:
            cont += 1
    return cont


PDB_PATH = path + "/data"
PATH_HOMOLWATDB = PDB_PATH + "/database/HomolWat_db.csv"

path_wat = PDB_PATH + "/PDB/REC_MOL_WAT/"
folder = glob.glob(path_wat + "*.pdb")

l_code, l_chain, l_naw, l_niw = [], [], [], []

for p2f in sorted(folder):
    file_name = p2f.split("/")
    code = file_name[-1][:-4]
    code2 = code.split("_")
    pdb = code2[0]
    chain = code2[1]
    rec = code2[2]
    # open file and save  CBs and wats
    f_in = open(p2f, "r")
    f_r = f_in.readlines()
    f_in.close()

    init_wats, CB_coords = [], []
    init_w_nums = []
    for line2 in f_r:
        res_name = line2[17:21].strip()
        at_type = line2[13:16].strip()
        res_num = line2[23:26].strip()
        if line2.startswith("ATOM"):
            if at_type == "CB":
                coxA = line2[30:38].strip()
                coyA = line2[38:46].strip()
                cozA = line2[47:54].strip()
                coord_cb = np.array([coxA, coyA, cozA], float)
                CB_coords.append(coord_cb)
        elif line2.startswith("HETATM"):
            if at_type == "O" and res_name == "HOH":
                init_wats.append(line2[:-1])
    if init_wats != []:
        prefilt_wats = []
        num_all_wats = len(init_wats)
        for wat in init_wats:
            Bfactor = float(wat[60:70].strip())
            if Bfactor < 45:
                prefilt_wats.append(wat)
        num_inn_wats = check_CircVar(CB_coords, prefilt_wats)
        
        l_code.append(code[:4])
        l_chain.append(code[5:6])
        l_naw.append(num_all_wats)
        l_niw.append(num_inn_wats)
    else:
        l_code.append(code[:4])
        l_chain.append(code[5:6])
        l_naw.append("0")
        l_niw.append("0")

# Create the table with the waters
data = {
    "PDB":l_code,
    "Chain":l_chain,
    "Water":l_naw,
    "Filtered waters": l_niw
}

waters_info = pd.DataFrame.from_dict(data)
waters_info.to_csv(PDB_PATH + '/database/waters_info.csv', index=False)

# Add the data into database: 
homolwat_db = pd.read_csv(PATH_HOMOLWATDB, sep=",")

homolwat_db_fil = pd.merge(homolwat_db, waters_info, left_on=['PDB accesion', 'Chain'], right_on=['PDB', 'Chain'], how='right')
homolwat_db_fil = homolwat_db_fil.drop(columns=['Waters'])
homolwat_db_fil = homolwat_db_fil.drop(columns=['PDB'])
print(homolwat_db_fil)
homolwat_db_fil.to_csv(PDB_PATH + '/database/HomolWat_fil_db.csv', index=False)
