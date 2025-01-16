#!/bin/bash -l
# Script to process protein structures, prepare models, superpose, and add waters.
# The script performs multiple steps and uses Python and PyMOL scripts for processing.

# Start time tracking
start=`date +%s`

# Help section
if [[ $1 == "-h" || $1 == "--help" ]]; then
    echo "Usage: $0 <protein_id> -s <state> -n <na_option>"
    echo
    echo "This script processes a protein structure through multiple steps:"
    echo "  1. Prepare folders and initial model"
    echo "  2. Superpose crystals to the model"
    echo "  3. Sort and group waters"
    echo "  4. Add internal water molecules or sodium ions"
    echo "  5. Merge the protein model with waters and renumber"
    echo "  6. Generate a PyMOL session with the final output"
    echo
    echo "Arguments:"
    echo "  <protein_id>   Identifier for the protein to process (e.g., 6oya)"
    echo "  -s <state>     Prioritize state: active, inactive, or intermediate (e.g., active)"
    echo "  -n <na_option> Specify 'with' or 'without' sodium ions (e.g., with)"
    echo
    echo "Output:"
    echo "  The final processed structure is saved as <protein_id>_HW.pdb in the output directory."
    echo
    echo "Example:"
    echo "  $0 6oya -s active -n with"
    exit 0
fi

# Check for mandatory inputs
if [[ -z $1 ]]; then
    echo "Error: No protein ID provided."
    echo "Use -h or --help for usage instructions."
    exit 1
fi

# Parse arguments
state=""
na_option=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--state) state="$2"; shift 2;;
        -n|--na) na_option="$2"; shift 2;;
        *) prot="$1"; shift;;
    esac
done

# Validate state input
if [[ -z $state ]]; then
    echo "Error: No state provided."
    echo "Use -h or --help for usage instructions."
    exit 1
fi

if [[ $state != "active" && $state != "inactive" && $state != "intermediate" ]]; then
    echo "Error: Invalid state '$state'. Must be 'active', 'inactive', or 'intermediate'."
    exit 1
fi

# Validate NA option
if [[ -z $na_option ]]; then
    echo "Error: No NA option provided."
    echo "Use -h or --help for usage instructions."
    exit 1
fi

if [[ $na_option != "with" && $na_option != "without" ]]; then
    echo "Error: Invalid NA option '$na_option'. Must be 'with' or 'without'."
    exit 1
fi

echo "Running Homolwat..." 
echo "- Protein ID: $prot"
echo "- State: $state"
echo "- NA Option: $na_option"

# MODIFY PATHS
# Set the working directory and define relevant paths
path=$(pwd)                            # Current working directory
path_rec=${path}/data/aln/             # Path to alignment data
path_wats=${path}/data/PDB/REC_WAT/    # Path to receptor water data
path_jobs=${path}/data/test/           # Path for job outputs
path_scripts=${path}/scripts/          # Path to the scripts directory
#path_pymol='/opt/soft/pymol-2/bin/pymol' # Uncomment and set if a custom PyMOL path is needed

# Input variables
pdb_folder=${prot:0:4}_HW              # Define folder name based on protein identifier
folder_out=$path_jobs$pdb_folder       # Output folder path

# Step 1: Prepare folders and MODEL (1/6)
module load BLAST+
echo 'Step 1) Preparing all..' 

python scripts/prep_all.py $prot $path $path_jobs $path_jobs $pdb_folder $path_scripts $state
echo 'Running Pymol part...'

# Step 2: Superpose crystals to the model (2/6)
echo 'Step 2) Superposing crystals to the model...' 
python scripts/prep_PDBs_mproc.py $pdb_folder $path_jobs $path_wats $path_scripts $pdb_folder

# Step 3: Sort and group waters (3/6)
echo 'Step 3) Sort and group waters...' 
python scripts/gene_wats_file_refined.py $pdb_folder $path_jobs $path_scripts

# Step 4: Add internal water molecules or sodium ion (4/6)
# Choose one of the following scripts:
echo 'Step 4) Add internal water molecules or sodium ion...' 
if [[ $na_option == "with" ]]; then
    python scripts/wat_adder_withNA.py $pdb_folder $path_jobs $path_scripts
else
    python scripts/wat_adder_noNA.py $pdb_folder $path_jobs $path_scripts
fi

# Step 5: Merge MODEL and waters, renumber waters, and prepare the output (5/6)
echo 'Step 5) Merging MODEL and waters, renumber waters, and prepare the output...' 
python scripts/merge_prot_waters.py ${folder_out}/HW/ ${pdb_folder}

# Step 6: Prepare a PyMOL session with the output (6/6)
echo 'Step 6) Preparing a PyMOL session with the output...' 
pymol -cqr scripts/gene_final_pse.py $pdb_folder $path_jobs $pdb_folder

# End time tracking and calculate execution time
end=`date +%s`
runtime=$((end-start))
echo "Execution time:" $runtime

# Cleanup: Delete temporary files
rm ${folder_out}/REC/*.pdb

# Final message: Show the location of the output
echo "Output is located at" ${folder_out} "as" ${pdb_folder}_HW.pdb