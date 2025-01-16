import glob
import random
import os
import sys

path = os.path.dirname(os.path.abspath("__file__"))
sys.path.insert(0, path)

Organism = {'BOVIN':'Bos taurus', 'MELGA': 'Meleagris gallopavo', 'HUMAN': 'Homo sapiens', 'TODPA': 'Todarodes pacificus',
            'MOUSE': 'Mus musculus', 'HCMVA': 'Cytomegalovirus', 'RAT': 'Rattus norvegicus', '9ARAC': 'Hasarius adansoni'}

def gener_table(archivo, table_id):
    path_f = path + "/data/REC_WAT/"
    table = ""
    pdbs, prot_names, clases, chains, organisms, resolutions, all_wats, in_wats, sodiums, ligands, dates, states = [], [], [], [], [], [], [], [], [], [], [], []
    botones, displais = [], []
    # table_in = open(folder_out + "/rec_sort_list.csv", "r")
    # table_r = table_in.readlines()
    for elems in sorted(archivo):
        # elems = line.split(" ")
        pdb, name_rec, clase, chain, organism, resol, all_wat, in_w, state, sodium, ligand, date = elems[0], elems[1], elems[2], elems[3], elems[4], elems[5], str(elems[6]), str(elems[7]),elems[11], elems[8], elems[9], elems[10]
        pdbs.append(pdb)
        prot_names.append(name_rec)
        clases.append(clase)
        chains.append(chain)
        organisms.append(organism)
        resolutions.append(resol)
        all_wats.append(all_wat)
        in_wats.append(in_w)
        sodiums.append(sodium)
        ligands.append(ligand)
        dates.append(date)
        states.append(state)
        display = "<input type='checkbox' class='display_box' value=" + pdb + "_" + chain + "_" + name_rec + " onclick=display()>"
        displais.append(display)
        boton =  "<a href=" +  path_f + pdb + "_" + chain + "_" + name_rec + ".pdb download><button>Download</button></a>"
        botones.append(boton)
    # table = ItemTable(items)
    table += '<table class="gpcr_table" id="' + table_id + '">\n'
    table += '<thead><tr><th>PDB</th><th>Receptor</th><th>Class</th><th>Chain</th><th>Organism</th><th>Resolution</th><th>All_wats</th><th>IN_wats</th><th>State</th><th>Sodium</th><th>Ligand</th><th>Date</th><th>Download</th><th>Display</th></tr></thead>\n'
    table += '<tbody>\n'
    for i in range(len(archivo)):
        table += '<tr><td><a href=https://www.rcsb.org/structure/' + pdbs[i] + '>' + pdbs[i] + '</a></td><td>' + prot_names[i] + '</td><td align=' + 'center' + '>' + clases[i] + '</td><td align=' + 'center' + '>' + chains[i] + '</td><td><i>' + organisms[i] + '</i></td><td align=' + 'center' + '>' + resolutions[i] + '</td><td align=' + 'center' + '>' + all_wats[i] + '</td><td align=' + 'center' + '>' + in_wats[i] + '</td><td align=' + 'center' + '>' + states[i] + '</td><td align=' + 'center' + '>' + sodiums[i] + '</td><td align=' + 'center' + '>' + ligands[i] + '</td><td>' + dates[i] + '</td><td>' + botones[i] + '</td><td align=' + 'center' + '>' + displais[i] + '</td></tr>\n'
        # <tr><td>ACM4</td><td>  23.1  </td><td>5DSG</td></tr>
    table += '</tbody>\n'
    table += '</table>\n'
    return table

path_l = path + "/scripts/"
lista = "gpcrListjul24.txt"
f_in = open(path_l + lista, "r")
f_r = f_in.readlines()

pdbs, prot_names, organisms, resolutions, dates = [], [], [], [], []
Info_PDB = {}
for line in f_r:
    elems = line.split(" ")
    pdb, name_rec, organism, clase, resol, date, state = elems[0], elems[1], Organism[elems[2]], elems[3], elems[6], elems[7], elems[8]
    Info_PDB[pdb] = [name_rec, organism, clase, resol, date, state]
    
path_all = path + "/REC_MOL_WAT/"
path_filt = path + "/REC_WAT/"

folder = glob.glob(path_all + "*.pdb")
final_list = []
names = []
for pdb_file in folder:
    all_wats, in_wats, sodium, ligand = 0, 0, 'No', ''
    path_parts = pdb_file.split('/')
    pdb_code = path_parts[-1]
    names.append(pdb_code[:-4])
    element = pdb_code.split("_")
    pdb = element[0]
    chain = element[1]
    rec_name = element[2][:-4]
    f_in_xtl = open(pdb_file, "r")
    f_r_xtl = f_in_xtl.readlines()
    for line in f_r_xtl:
        if line.startswith('HETATM'):
            molec = line[17:20].strip()
            if molec == "HOH":
                all_wats += 1
    f_in_filt = open(path_filt + pdb_code[:-4] + ".pdb", "r")
    f_r_filt = f_in_filt.readlines()
    tmp_lig = []
    for line2 in f_r_filt:
        if line2.startswith('HETATM'):
            molec = line2[17:20].strip()
            if molec == "NA":
                sodium = "Yes"
            if molec == "HOH":
                in_wats += 1
            if molec != 'NA' and molec != 'HOH':   
                if molec not in tmp_lig:
                    # print molec
                    tmp_lig.append(molec)
    # print pdb, tmp_lig
    if len(tmp_lig) > 0:
        ligand = "Yes"
    else:
        ligand = "No"
        # if len(tmp_lig) == 1:
        #     ligand = tmp_lig[0]
        # else:
        #     ligand = " ".join(tmp_lig)
    info_out = [pdb, rec_name, Info_PDB[pdb][2], chain, Info_PDB[pdb][1], Info_PDB[pdb][3], str(all_wats), str(in_wats), sodium, ligand, Info_PDB[pdb][4], Info_PDB[pdb][5]]
    # print " ".join(info_out)
    final_list.append(info_out)
# 
# genera tabla gpcr

f_out_gpcr = open(path + "/data/web/table_gpcr.html", "w")
f_out_gpcr.write('<script src="static/datatable_gpcr.js"></script>\n')
f_out_gpcr.write(gener_table(final_list, 'gpcr_table'))
f_out_gpcr.close()

# genera tabla view

f_out_view = open(path + "/data/web/table_view.html", "w")
f_out_view.write('<script src="static/datatable_view.js"></script>\n')
f_out_view.write(gener_table(final_list, 'view_table'))
f_out_view.close()

#print gener_table(final_list)
## generar diccionarios de colores
d_pdb_color = {}
for d in names:
    r = lambda: random.randint(0, 255)
    d_pdb_color[d] = "#%02X%02X%02X" % (r(), r(), r())

# 3. Copias y pegas el output en static/ngl_data.html por la variable var d_pdb_color = {...}
colors = "var d_pdb_color={"
# print(colors)
for k, v in d_pdb_color.items():
    s = "'" + k + "': '" + v + "'"
    if colors[-1] == "{":
        colors += s
    else:
        colors += "," + s
colors_dict = colors + '}'
f_out_col = open(path + "/data/web/colors.js", "w")
f_out_col.write(colors_dict)

