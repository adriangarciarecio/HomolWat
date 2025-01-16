#!/usr/bin/env python3

import pymol
import os
import glob
import sys
#import a_pdb2ball
from clize import run
import wget
import urllib
import requests

import pandas as pd

path = os.path.dirname(os.path.abspath("__file__"))
sys.path.insert(0, path)

allres = "ALA+ARG+ASP+GLU+THR+TYR+HIS+LEU+ILE+VAL+PRO+GLY+TRP+LYS+PHE+ASN+GLN+MET+CYS+SER+ACE+NH2+DI7+DI8"
wat_ions = "HOH+NA+CL+ZN+SO4"
dirs = ["FULL", "CHAIN", "REC_MOL_WAT"]

#PREPARE THE CHAINS, UNIPROT ACC, EXP 
def get_pdb_chains():
    """Gets the list of PDBs and their chains from PDBLIB"""
    all_info_rec = pd.read_csv(PATH_HOMOLWATDB, sep=",")
    pdb_list = all_info_rec["PDB accesion"].unique().tolist()
    d_chains = {}
    for pdb in pdb_list: 
        f_chains = all_info_rec[all_info_rec["PDB accesion"] == pdb]
        d_chains[pdb] = f_chains["Chain"].to_list()
    #d_chains = dict(list(zip(all_info_rec["PDB accesion"], all_info_rec["Chain"])))
    d_uniprot = dict(list(zip(all_info_rec["PDB accesion"], all_info_rec["Uniprot entry"])))
    d_expmet = dict(list(zip(all_info_rec["PDB accesion"], all_info_rec["Method"])))
    d_segment = dict(list(zip(all_info_rec["PDB accesion"], all_info_rec["Segments Pdb"])))
    return pdb_list, d_chains, d_uniprot, d_expmet, d_segment

def prepare_FULL(pdb_list):
    """Prepares the direcory full.
    PDB_LIST is the list of PDBs in PDBLIB;
    pdb_pepts contains the chain mappings;
    uniprot contains the ids from PDBLIB"""
    # import shutil
    import os

    new_pdbs = []  # xxxx_A
    os.chdir(f"{PDB_PATH}/PDB/FULL")
    for pdb in pdb_list:  # these are xxxx(.pdb)
        if str(pdb) == "nan":
            continue
        print(pdb, "in pdb_list")
        # DO SOMETHING ONLY IF PDB IS NEW!!!
        if (
            (not os.path.exists(pdb.lower() + "_final.pdb"))
            and (not os.path.exists(pdb.lower() + "_updated.cif"))
            and (not os.path.exists(pdb.lower() + ".cif"))
        ):
            found_pdb = False
            if d_expmet[pdb] == "X-ray":
                try:
                    url = (
                        f"https://pdb-redo.eu/db/{pdb.lower()}/{pdb.lower()}_final.pdb"
                    )
                    wget.download(url)
                    found_pdb = True
                    print(f"\nDownloading PDB-REDO structure {pdb}\n")
                except urllib.error.HTTPError:
                    pass
            if not found_pdb:
                try:
                    url = f"http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb.lower()}_updated.cif"
                    wget.download(url)
                    found_pdb = True
                    print(f"\nDownloading updated CIF PDBe structure {pdb}\n")
                except urllib.error.HTTPError:
                    try:
                        url = f"http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb.lower()}.cif"
                        wget.download(url)
                        found_pdb = True
                        print(f"\nDownloading updated CIF PDBe structure {pdb}\n")
                        # http://www.ebi.ac.uk/pdbe/entry-files/download/pdb2rh1.ent
                        # prevents strange errors when a PDB cannot be retrieved (it will the next run)
                        #
                    except urllib.error.HTTPError:
                        print("Not found!")
                        sys.exit()
            if found_pdb:
                new_pdbs.append(pdb)
        else:
            print("Already processed. Remove file from ./FULL to regenerate it")

    # ends in FULL
    os.chdir("..")
    print(20 * "-" + "\n")
    print("        --- FINISHED --- directory FULL ---")
    return new_pdbs

def prepare_CHAIN(new_pdbs, d_chains, d_uniprot):
    import glob

    print(new_pdbs)
    os.chdir(f"{PDB_PATH}/PDB/FULL")
    pdbs_in_cif = [pdb[:4] for pdb in glob.glob("*.cif")]
    print("pdbs_in_cif", pdbs_in_cif)
    pdbs_from_pdbredo = [pdb[:4] for pdb in glob.glob("*_final.pdb")]
    print("pdbs_from_pdbredo", pdbs_from_pdbredo)

    os.chdir(f"{PDB_PATH}/PDB")
    for pdb in new_pdbs:  # these are 1f88_updated.cif
        print(pdb)
        pdb = pdb[0:4].lower()
        # print ' Chains', d_chains[pdb]
        pymol.cmd.do("cd " + "./FULL")
        print("<><>" + 20 * "-" + "\n" + pdb)
        if pdb.lower() in pdbs_in_cif:
            pymol.cmd.do("load " + pdb + "_updated.cif, " + pdb)
            pymol.cmd.do(f"save {pdb}.pdb, {pdb}")
            pymol.cmd.do("delete " + pdb)
            # os.remove(pdb + ".cif")
            pymol.cmd.do("load " + pdb + ".pdb")
        elif pdb.lower() in pdbs_from_pdbredo:
            pymol.cmd.do("load " + pdb + "_final.pdb, " + pdb)
            pymol.cmd.do(f"save {pdb}.pdb, {pdb}")
            pymol.cmd.do("delete " + pdb)
            # os.remove(pdb + "_final.pdb")
            pymol.cmd.do("load " + pdb + ".pdb")
        else:
            print("No PDB to load")
            sys.exit()
        pymol.cmd.do("cd " + "../CHAIN")
        #l_ch = d_chains[pdb.upper()].split(",")
        l_ch = d_chains[pdb.upper()]
        for ch in l_ch:
            # print(ch)
            # Break the PDB in CHAINs
            # This is for peptides, which have their own chain!!
            # with_pept = False
            # for pdb_pept in pdb_pepts:
            #     if pdb in pdb_pept[0]:
            #         with_pept = True
            #         print(pdb_pept[1])  # ????
            #         if ch in pdb_pept[1]:  # avoid unknown chains as keys
            #             cmd_txt = "select sele_ch, chain {}+{}".format(
            #                 pdb_pept[1][ch], ch
            #             )
            #             pymol.cmd.do(cmd_txt)
            #             break
            print("----->", pdb)
            # THE NEXT FOUR LINES WERE INTRODUCED RECENTLY AND HAVE NOT BEEN TESTED
            # pymol.cmd.do("select prot, resn " + allres + " and " + pdb.lower())
            # pymol.cmd.do("select close, (br. prot around 5) and " + pdb.lower())
            # pymol.cmd.do("select far, not (byres (close or prot) around 5) and not (close or prot) and " + pdb.lower())
            # pymol.cmd.do("remove far")
            # if not with_pept:  # general case
            pymol.cmd.do("select sele_ch, chain " + ch)
            # print pdb, ch, d_uniprot[pdb]  # for testing
            pymol.cmd.do(
                "save " + pdb.upper() + "_" + ch + "_" + d_uniprot[pdb.upper()] + ".pdb, sele_ch"
            )
            if pdb.upper() == "4EIY" and ch == "A":
                pymol.cmd.do(
                "save ref_4EIY_A.pdb, sele_ch"
                )

        pymol.cmd.do("delete *")
        pymol.cmd.do("cd " + "..")

    # ends in FULL
    print(20 * "-" + "\n")
    print("        --- FINISHED --- directory CHAINS ---")

# Check
def prepare_REC(new_pdbs, d_segments):
    """Actually prepares  REC_MOL_WAT....
    Prepare PDBs with/without WAT, MOL all supperposed to "refpdb"
    "new_pdbs" is a list"""
    l_error_wat = []
    print(os.getcwd())
    os.chdir(f"{PDB_PATH}/PDB/CHAIN")
    pdb_c = glob.glob("????_*_*.pdb")
    #pdb_c = ['2RH1_A_ADRB2.pdb', '4DKL_A_OPRM.pdb']  # for testing
    rms_a, rms_s, rms_pdb = [], [], []
    refpdb = "ref_4EIY_A"  # this is the reference structure for superposition
    pymol.cmd.do("load " + refpdb + ".pdb")
    os.chdir(f"{PDB_PATH}/PDB")

    # Loop stats at dir MAIN
    #print("new_pdbs", new_pdbs)
    d_waters = {}
    for pdbname in pdb_c:  # these are xxxx_A.pdb
        shrt_pdbn = pdbname[:-4]  # strip extension
        # if 1 ==1:  # for testing; replace if
        #if shrt_pdbn[:4].upper() in new_pdbs:  # for comparison take only pdbid
        print("#############", shrt_pdbn, pdbname)
        pymol.cmd.do("cd " + "./CHAIN")
        pymol.cmd.do("load " + pdbname)
        print(pdbname)
        # save 1 single chain without lysozyme, BRIL ....
        # EXCEPTIONs: BRIL/T4L at the N-terminus ...
        chain = pdbname[5]
        segment = d_segments[pdbname[0:4]]
        segment = segment.replace(",", "+")
        print(segment)
        try:
            #select not resi 1-352 and (segid A and 4DKL_A_OPRM_MOUSE)
            pymol.cmd.do("select solvent and "
                + shrt_pdbn)
            wat = pymol.cmd.count_atoms("sele")
            print(wat)
            d_waters[pdbname[0:4]] = int(wat)
            print(d_waters)
            pymol.cmd.do(
                "select not resi " 
                + segment
                + " and not hetatm and "
                + shrt_pdbn 
            )  # before 900-5000
            pymol.cmd.do("remove sele")
            # remove water + ion that are very far
            pymol.cmd.do("select prot, resn " + allres + " and " + shrt_pdbn)
            pymol.cmd.do("select close, (br. prot around 5) and " + shrt_pdbn)
            pymol.cmd.do(
                "select far, not (byres (close or prot) around 5) and not (close or prot) and "
                + shrt_pdbn
            )
            pymol.cmd.do("remove far")
            # Superposition of structures
            rms_data_super = pymol.cmd.super(
                shrt_pdbn + "////CA" + " and chain " + chain, refpdb + "////CA"
            )
            rms_data_align = pymol.cmd.align(shrt_pdbn + " and chain " + chain, refpdb)
            rms_super = rms_data_super[0]
            rms_align = rms_data_align[0]
            print("@@@@@ RMS super, align", rms_data_super[0], rms_data_align[0])
            print(shrt_pdbn[:4])
            # With exceptions to SUPER giving larger RMS than ALI but being better
            if rms_align > rms_super or shrt_pdbn[:4] == "6N51":
                rms_data_super = pymol.cmd.super(
                    shrt_pdbn + "////CA", refpdb + "////CA"
                )
            rms_a.append(rms_align)
            rms_s.append(rms_super)
            rms_pdb.append(shrt_pdbn)
            pymol.cmd.do("cd " + "../REC_MOL_WAT")
            pymol.cmd.do("save " + pdbname + ", " + shrt_pdbn)
        except:
            print("Error: " + shrt_pdbn)
            pymol.cmd.do("cd " + "../REC_MOL_WAT")
        # end loop
        l_error_wat.append(shrt_pdbn)
        pymol.cmd.do("delete " + shrt_pdbn)
        pymol.cmd.do("cd ..")
    print(d_waters)
    print("        --- FINISHED --- directories REC_MOL_WAT, REC_WAT, REC  ---")

    # write waters to db
    all_info_rec = pd.read_csv(PATH_HOMOLWATDB, sep=",")
    all_info_rec_clean = all_info_rec.copy()
    l_waters = []
    for index, row in all_info_rec.iterrows(): 
        try:
            pdb_n = row["PDB accesion"] 
            l_waters.append(d_waters[pdb_n.upper()])
        except Exception as e:
            print(f"Error at index {index}: {e}")
            all_info_rec_clean.drop(index, inplace=True)
    all_info_rec_clean["Waters"] = l_waters
    
    # write fasta to db 
    l_fasta = []
    l_uni = all_info_rec_clean["Uniprot accesion"].to_list()
    for uni in  l_uni:
        url_uni = f'https://www.uniprot.org/uniprot/{uni}.fasta'
        response = requests.get(url_uni)
        data_uni = response.text
        l_fas = data_uni.split("\n")
        fas = "".join(l_fas[1:])
        l_fasta.append(fas)
    all_info_rec_clean["Fasta"] = l_fasta
    all_info_rec_clean.to_csv(PATH_HOMOLWATDB, index=False)

    # write rmsds
    file_rms = open("rmsd.txt", "w")
    rms_sort = list(zip(rms_pdb, rms_s, rms_a))
    rms_sort.sort()
    # rms_sort.reverse()
    print("\n", rms_sort, "\n")
    for rms in rms_sort:
        out = "%s %6.3f %6.3f\n" % (rms[0], rms[1], rms[2])
        file_rms.write(out)
    file_rms.close()
    return l_error_wat

def prepare_REC_MOL(new_pdbs):
    # Prepare PDBs REC_MOL
    os.chdir("./data/PDB/REC_MOL_WAT")
    # pymol.cmd.do("cd " + '../REC_MOL_WAT')
    pdb_c = glob.glob("????_*_*.pdb")
    for pdbname in pdb_c:
        # print pdbname
        shrt_pdbn = pdbname[:-4]
        if shrt_pdbn[:4] in new_pdbs:
            # print shrt_pdbn
            pymol.cmd.do("load " + pdbname)
            pymol.cmd.do("cd " + "../REC_MOL")
            pymol.cmd.do("select resn " + wat_ions + " and " + shrt_pdbn)
            pymol.cmd.do("remove sele")
            pymol.cmd.do("save " + pdbname + ", " + shrt_pdbn)
            pymol.cmd.do("delete " + shrt_pdbn)
            pymol.cmd.do("cd " + "../REC_MOL_WAT")
    print("        --- FINISHED --- directories REC_MOL  ---")

################## MAIN PROGRAM #############################################################

def main():
    """Download and process PDBs associated to "lmcdb pdb"
    It will only process new PDBs.
    If you want to regenerate old PDBs remove them from the directory FULL.
    problem

    # For new proteins, check (in the directory FULL) which are the chains that map one receptor.
    # Ex: receptor + peptide, receptor + nanobody ...
    # In case there is more than "receptor" , edit "pdb_mappings.txt", remove the new PDBs and rerun this
    :param machine: where to run the sctip (uab or arnau)
    """
    global host, PDB_PATH, PATH_HOMOLWATDB, pdb_pepts
    global pdb_list, d_chains, d_uniprot, d_expmet

    #DEFINE PATHS
    PDB_PATH = path + "/data"
    PATH_HOMOLWATDB = PDB_PATH + "/database/HomolWat_db.csv"

    print(os.getcwd())
    if not os.path.exists(f"{PDB_PATH}/PDB"):
        os.mkdir(f"{PDB_PATH}/PDB")
    os.chdir(f"{PDB_PATH}/PDB")

    if not os.path.exists("CHAIN"):
        for dir in dirs:
            try:
                os.stat(dir)
            except FileNotFoundError:
                os.mkdir(dir)

    pdb_list, d_chains, d_uniprot, d_expmet, d_segments = get_pdb_chains()
    # for testing
    # pdb_list = ['1GZM', '2RH1', '4DKL']
    # d_chains = {'1GZM': ['A', 'B'], '2RH1': ['A'], '4DKL': ['A'],}
    # d_uniprot = {'1GZM': 'OPSD', '2RH1': 'ADRB2','4DKL': 'OPRM',}

    #pdb_pepts = a_pdb2ball.read_pdb_mapping()
    new_pdbs = prepare_FULL(pdb_list)  # returns the list of only "new" pdb
    # pymol.finish_launching()  # maybe ['pymol', '-q'])
    prepare_CHAIN(new_pdbs, d_chains, d_uniprot)  # CHANGE new_pdbs, pdb_list
    l_error_wat = prepare_REC(new_pdbs, d_segments)  # CHANGE new_pdbs, pdb_list
    # prepare_REC_MOL(new_pdbs)

    pymol.cmd.quit()
    print(l_error_wat)
    
if __name__ == "__main__":
    run(main)
