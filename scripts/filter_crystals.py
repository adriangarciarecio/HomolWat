"""
from cristals with internal waters
clean and keep just:
    - orthosteric ligand
    - internal water and sodium
"""

import math
import numpy as np
import glob
from EM_functions import *
import os
import sys

path = os.path.dirname(os.path.abspath("__file__"))
sys.path.insert(0, path)
print(path)
path_in = path + '/data/PDB/REC_MOL_WAT/'
path_out = path + "/data/PDB/REC_WAT/"

folder = glob.glob(path_in + "*.pdb")

for p2f in sorted(folder):
    print(p2f)
    file_name = p2f.split("/")
    code = file_name[-1][:-4]
    # print(code)
    code2 = code.split("_")
    pdb = code2[0]
    chain = code2[1]
    rec = code2[2]
    print (pdb, chain, rec)
    ### open file and save  CBs and wats
    f_in = open(p2f, "r")
    f_r = f_in.readlines()
    f_in.close()

    init_wats, CB_coords = [], []
    init_w_nums, protein, ions, ligands = [], [], [], []
    for line2 in f_r:
        res_name = line2[17:21].strip()
        at_type = line2[13:16].strip()
        res_num = line2[23:26].strip()
        if line2.startswith("ATOM"):
            protein.append(line2)
            if at_type == "CB":
                coxA = line2[30:38].strip()
                coyA = line2[38:46].strip()
                cozA = line2[47:54].strip()
                coord_cb = np.array([coxA, coyA, cozA], float)
                CB_coords.append(coord_cb)
        elif line2.startswith("HETATM"):
            if res_name == "HOH":
                init_wats.append(line2)
            elif res_name == "NA":
                ions.append(line2)
            else:
                ligands.append(line2)
    filt_wats = []
    if init_wats != []:
        num_inn_wats = 0
        prefilt_wats = []
        num_all_wats = len(init_wats)
        for wat in init_wats:
            coxA = wat[30:38].strip()
            coyA = wat[38:46].strip()
            cozA = wat[47:54].strip()
            coord_w = np.array([coxA, coyA, cozA], float)
            # coords.append(coord_i)
            Bfactor = float(wat[60:70].strip())
            if Bfactor < 45:
                # prefilt_wats.append(wat)
                if check_CircVar(CB_coords, coord_w, 10) > 0.6:
                    num_inn_wats += 1
                    filt_wats.append(wat)
        print (code, num_all_wats, num_inn_wats)
    int_mol, ligand_lines = [], []
    if ligands != []:
    #vamos a guardar solo el ligando interno
        # 1o hacemos lista de los ligandos del pdb
        molecs = []
        for l in ligands:
            molec = l[17:20].strip()
            num = l[22:27].strip()
            info = [molec, num]
            if info not in molecs:
                molecs.append(info)
        # 2o sacamos el COM de cada uno
        for i, mol in enumerate(molecs):
            coords = []
            #print i, mol
            for l2 in ligands:
                molec = l2[17:20].strip()
                num = l2[22:27].strip()
                if molec == mol[0] and num == mol[1]:
                    #print l2
                    coxA = l2[30:38].strip()
                    coyA = l2[38:46].strip()
                    cozA = l2[47:54].strip()
                    coord_i = np.array([coxA, coyA, cozA], float)
                    coords.append(coord_i)
            com_mol = center_of_molec(coords)
            molecs[i].append(com_mol)
        # 3o sacamos el CV de la molecula
        for j, mol in enumerate(molecs):
            cv = check_CircVar(CB_coords, mol[2], 20)
            molecs[j].append(cv)
        #for x in molecs:
            #print x
    # 4o ordenamos la lista de molecs por cv y nos quedamos el ultimo (max CV)
        sorted_mols = sorted(molecs, key=lambda cv: cv[3])
        #print "ordenados"
        int_mol_num = [sorted_mols[-1][0], sorted_mols[-1][1]]
        #print int_mol_num
        # 5o se guardan las lineas del ligando en una lista
        # if len(molecs) > 1:
        #     f_out_other_ligands = open(path_HW + "/other_ligs.pdb", "w")
        for line in ligands:
            molec = line[17:20].strip()
            num = line[22:27].strip()
            if molec == int_mol_num[0] and num == int_mol_num[1]:
                #print line
                ligand_lines.append(line)



    f_out = open(path_out + code + ".pdb", "w")
    for l1 in protein:
        f_out.write(l1)
    for l4 in ligand_lines:
        f_out.write(l4)
    if ions != []:
        for l2 in ions:
            f_out.write(l2)
    for l3 in filt_wats:
        f_out.write(l3)
        
