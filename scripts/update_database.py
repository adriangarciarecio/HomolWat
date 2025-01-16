# MAIN IMPORT
import requests
import pandas as pd
import json
import os 
from time import sleep
import sys
import subprocess
from datetime import datetime
import re 

path = os.path.dirname(os.path.abspath("__file__"))
sys.path.append(os.path.dirname(path))

# IMPORTS
from SECRETS import *
from pymol import cmd
#from pypdb import *


# DATABASE IMPORTS
from sqlalchemy import create_engine


def get_opm_coords(pdb):
    '''Get the coordinates to determine the differents zones on the structure of a protein.'''

    # try:
    response = requests.get(f"https://lomize-group-opm.herokuapp.com/assemblies/pdb/{pdb}")
    print(response)
    opm_data = response.json()
    opm_id = opm_data["id"]
    pdb_id = opm_data["pdbid"]
    seq = opm_data["aaseq"]
    tm_seg = opm_data["transmembrane_segment"]
    tm_count = opm_data["transmembrane_alpha_helix_count"]
    data = {
        "opm_id": opm_id,
        "pdb": pdb_id,
        "seq":seq,
        "tm_seg":tm_seg,
        "tm_count":tm_count,
    }
    return data
    # except Exception as e:
    #     print(f"- Error {e} on the request of the PDB {pdb}")

def update_opm_info():
    # Run the connection engine to the database
    print("- Engine on the connection to MemProtDb database")
    engine = create_engine(f'postgresql+psycopg2://{DB_USER}:{DB_PASSWORD}@localhost/{DB_NAME}')
    
    # Get list of pdbs
    sql_query = f"SELECT pdb, uniprot_accesion FROM gpcrdb_info WHERE resolution < 90" #and type='X-ray diffraction';
    data = pd.read_sql(sql_query, engine)
    
    pdbids = data["pdb_code"].values
    uniprot = data["protein"].values
    print(uniprot)
    # pdbids = ["4DKL"]
    
    print(f"- Getting information from OPM")
    opm_dict = {
        "opm_id":[],
        "pdb":[],
        "seq":[],
        "tm_seg":[],
        "tm_count":[],
        
    }
    for i, uni in enumerate(uniprot):
        print(f"- PDB: {pdbids[i]}")
        opm_pdb_data = get_opm_coords(uni.lower())
        print(opm_pdb_data)
        opm_dict["opm_id"].append(opm_pdb_data["opm_id"])
        opm_dict["uniprot"].append(uni)
        opm_dict["pdb"].append(pdbids[i])
        opm_dict["seq"].append(opm_pdb_data["seq"])
        opm_dict["tm_seg"].append(opm_pdb_data["tm_seg"])
        opm_dict["tm_count"].append(opm_pdb_data["tm_count"])

    # Store the info in the database
    print("- Update opm_info table")
    opm_info = pd.DataFrame(opm_dict)
    print(opm_info)
    opm_info.to_sql(name='opm_info', if_exists='replace', con=engine)

def update_gpcrdb_info():
    
    # Get data from GPCRdb 
    print("- Getting information from GPCRdb")
    strucgpcrdb=requests.get('https://gpcrdb.org/services/structure/').json()
    gpcrdb_info = pd.DataFrame.from_dict(pd.json_normalize(strucgpcrdb), orient='columns')
    gpcrdb_info = gpcrdb_info.sort_values(by=['pdb_code'])
    gpcrdb_info["ligands"] = gpcrdb_info["ligands"].apply(json.dumps)
    print(gpcrdb_info)
    
    # Run the connection engine to the database
    print("- Engine on the connection to MemProtDb database")
    engine = create_engine(f'postgresql+psycopg2://{DB_USER}:{DB_PASSWORD}@localhost/{DB_NAME}')
    
    # Store the info in the database
    print("- Update gpcrdb_info table")
    gpcrdb_info.to_sql(name='gpcrdb_info', if_exists='replace', con=engine)

def update_pdb_info():  
    
    # PFAM GPCRs
    pfamdict = {
    "classA": "PF00001",
    "classB": "PF00002",
    "classC": "PF00003",
    "frizzled": "PF01534",}
    
    # Run the connection engine to the database
    print("- Engine on the connection to MemProtDb database")
    engine = create_engine(f'postgresql+psycopg2://{DB_USER}:{DB_PASSWORD}@localhost/{DB_NAME}')
    
    # Get list of pdbs
    sql_query = f"SELECT pdb_code, preferred_chain, protein FROM gpcrdb_info WHERE resolution < 90" #and type='X-ray diffraction';
    data = pd.read_sql(sql_query, engine)
    
    pdbids = data["pdb_code"].values
    
    # test
    # pdbids = ["8C9W", "8GNG"]
        
    all_info, l_seg_pdb, l_error, l_seg = [], [], [], []
    clas = ""
    j = 1
    for pdb in pdbids:
        print(f"> {j *100 / len(pdbids)} %") 
        print(f"- Getting information for {pdb}")
        url_pdb=f"https://data.rcsb.org/rest/v1/core/entry/{pdb.lower()}"
        response = requests.get(url_pdb)
        data = json.loads(response.text)
        try:
            res = data["rcsb_entry_info"]["diffrn_resolution_high"]["value"]
        except:
            try:
                res = data["rcsb_entry_info"]["resolution_combined"][0]
            except:
                l_error.append(pdb) 
                continue
        try: # When data not found 
            meth = data["rcsb_entry_info"]["experimental_method"]
        except:
            l_error.append(pdb) 
            continue
        ids = data["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
        for id in ids:
            url = f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb.lower()}/{id}'
            response = requests.get(url)
            data = json.loads(response.text)
            #chain = data["entity_poly"]["pdbx_strand_id"]# la chain que correspon al polymer entity  
            chain_pdb = data["rcsb_polymer_entity_container_identifiers"]["asym_ids"]
            chain_auth = data["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]
            for i, ch in enumerate(chain_pdb):
                mapping = f'https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb.lower()}/{ch}' 
                response = requests.get(mapping)
                data_map = json.loads(response.text)
                auth_map = data_map["rcsb_polymer_entity_instance_container_identifiers"]["auth_to_entity_poly_seq_mapping"]# ["852","853","854","855","856","
                try:
                    l_rcsb = data["rcsb_polymer_entity_align"]
                    for rcsb in l_rcsb: 
                        # uniprot
                        acc=rcsb['reference_database_accession']
                        url_uni = f'https://rest.uniprot.org/uniprotkb/{acc}.txt'
                        response = requests.get(url_uni)
                        data_uni = response.text
                        # Get the Uniprot entry
                        search = re.search('ID   ()\w+', data_uni)    
                        entry = str(search.group()).replace("ID", "").replace(" ", "").lower()
                        # extreiem el pfam 
                        search = re.search("Pfam; \w+", data_uni)    
                        pfam = str(search.group()).replace("Pfam;", "").replace(" ", "")
                        if pfam in pfamdict.values():
                            clas = "receptor"
                        else:
                            clas = "fusion protein"
                        #Segment region
                        l_rscb_fea = data["rcsb_polymer_entity_feature"] 
                        for rscb_fea in l_rscb_fea:
                            # el receptor comença un segment a ref_beg_seq_id 
                            l_reg = rcsb['aligned_regions']
                            for reg in l_reg:
                                res_ini_pdb = reg['entity_beg_seq_id'] - 1  # 1 --: 0 
                                res_end_pdb = res_ini_pdb + reg['length'] - 1 
                                res_ini_auth = auth_map[res_ini_pdb] 
                                res_end_auth = auth_map[res_end_pdb] 
                                l_seg.append(f"{res_ini_pdb}-{res_end_pdb}")
                                l_seg_pdb.append(f"{res_ini_auth}-{res_end_auth}")
                            if not [pdb,res,meth,id,acc,entry,pfam,chain_auth[i][0],",".join(l_seg_pdb), clas] in all_info:
                                if clas != "":
                                    all_info.append([pdb,res,meth,id,acc,entry,pfam,chain_auth[i][0],",".join(l_seg_pdb), clas])
                                elif clas == "" and pfam == "":
                                    all_info.append([pdb,res,meth,id,acc,entry,pfam,chain_auth[i][0],",".join(l_seg_pdb), clas])
                            l_info, l_seg_pdb, l_seg = [], [], []
                            clas = ""
                except Exception as e:
                    if not [pdb,res,meth,id,"","","",chain_auth[i][0],"","", ""] in all_info:
                        all_info.append([pdb,res,meth,id,"","","",chain_auth[i][0],"","", ""])
                    l_info, l_seg_pdb, l_seg = [], [], []
                    clas = ""
                    if pdb not in l_error:
                        l_error.append(pdb) 
        j += 1
        
    print(l_error)
    gpcr_pdb = pd.DataFrame(all_info, columns=['pdb', 'resolution', 'method', 'id', 'uniprot_accesion', 'uniprot_entry', 'pfam', 'chain', 'segments', 'type', 'what'])
    gpcr_pdb = gpcr_pdb.drop(['what'],axis = 1)
    print(gpcr_pdb)
    
    # Store the info in the database
    print("- Update gpcrdb_info table")
    gpcr_pdb.to_sql(name='gpcr_pdb', if_exists='replace', con=engine)
    
# MAIN
def backup_database(username, database, save_path, backup_file, host='localhost', port='5432'):
    # Construct the pg_dump command
    command = [
        'pg_dump',
        '-U', username,
        '-h', host,
        '-p', port,
        database
    ]

    # Redirect output to the backup file
    with open(f"{save_path}/{backup_file}", 'w') as f:
        try:
            subprocess.run(command, stdout=f, check=True, env={"PGPASSWORD": DB_PASSWORD})
            print(f"Backup of database '{database}' completed successfully and saved to '{save_path}/{backup_file}'.")
        except subprocess.CalledProcessError as e:
            print(f"Error occurred during backup: {e}")

if __name__ == "__main__":
     # Get the current timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    SAVE_PATH = sys.argv[1]
    BACKUP_FILE = f'backup_{DB_NAME}_{timestamp}.sql'

    # Execute the backup
    update_gpcrdb_info()
    update_pdb_info()
    # update_opm_info()
    backup_database(DB_USER, DB_NAME, SAVE_PATH, BACKUP_FILE)
