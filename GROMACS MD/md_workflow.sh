#!/bin/bash

################################################################################
# Molecular Dynamics Simulation Workflow (Mac Version - Auto Selection)
# Based on CHARMM36 Force Field and GROMACS
# Automatically selects force fields and groups
################################################################################

set -e  # Exit on any error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log() {
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

################################################################################
# Initialize Conda (Mac/Linux compatible)
################################################################################
init_conda() {
    log "Initializing conda..."
    
    # Determine shell type and initialize conda accordingly
    if [[ $SHELL == *"zsh"* ]]; then
        eval "$(conda shell.zsh hook)" 2>/dev/null || true
    else
        eval "$(conda shell.bash hook)" 2>/dev/null || true
    fi
}

################################################################################
# PREREQUISITES CHECK
################################################################################
check_prerequisites() {
    log "=== CHECKING PREREQUISITES ==="
    
    # Check for CHARMM36
    if [[ ! -d "charmm36-jul2022.ff" ]]; then
        error "charmm36-jul2022.ff directory not found"
        echo "Download from: http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs"
        exit 1
    fi
    success "CHARMM36 force field found"
    
    # Check for conda
    if ! command -v conda &> /dev/null; then
        error "conda not found. Please install miniconda3"
        exit 1
    fi
    success "Miniconda3 found"
    
    # Initialize conda for this shell
    init_conda
    
    # Check md_core environment exists
    if ! conda env list | grep -q "md_core"; then
        error "md_core environment not found"
        echo "Create it with: conda create -n md_core -c conda-forge gromacs pdbfixer numpy"
        exit 1
    fi
    success "md_core environment found"
    
    # Activate and check for GROMACS
    conda activate md_core 2>/dev/null
    if ! command -v gmx &> /dev/null; then
        error "GROMACS not found in md_core"
        echo "Install with: conda activate md_core && conda install -c conda-forge gromacs"
        exit 1
    fi
    success "GROMACS found"
    
    log "All prerequisites verified!"
}

################################################################################
# PDB PREPARATION
################################################################################
prepare_pdb() {
    log "=== PDB PREPARATION ==="
    
    local pdb="${1:-FtsZ-Dacomitinib128.pdb}"
    
    if [[ ! -f "$pdb" ]]; then
        error "PDB file not found: $pdb"
        return 1
    fi
    
    log "Please manually change all 'UNL' to 'LIG' in: $pdb"
    read -p "Press Enter when done, or 'q' to skip: " response
    [[ "$response" == "q" ]] && return 1
    
    success "PDB file ready"
}

################################################################################
# SEPARATE PDB
################################################################################
separate_pdb() {
    log "=== SEPARATING PDB ==="
    
    local pdb="${1:-FtsZ-Dacomitinib128.pdb}"
    
    init_conda
    conda activate md_core 2>/dev/null
    
    log "Extracting protein..."
    grep -E "ATOM|HETATM" "$pdb" | grep -v "LIG" > protein_GDP.pdb
    
    log "Extracting ligand..."
    grep "LIG" "$pdb" > LIG.pdb
    
    log "Adding crystal structure info..."
    {
        echo "CRYST1   67.715   54.650   84.496  90.00 106.75  90.00 C 1 2 1       1"
        awk '
        /^TER/ { next } # remove old TER
        /^END/ { next } # remove old END

        /^HETATM/ && !ter_added {
            print "TER"
            ter_added=1
        }

        { print }

        END {
            if (!ter_added) {
                print "TER" # fallback: ensure TER exists
            }
        
        }
        ' protein_GDP.pdb 
        echo "END"
    } > protein_GDP.tmp && mv protein_GDP.tmp protein_GDP.pdb
    
    log "Running pdbfixer..."
    pdbfixer protein_GDP.pdb --add-atoms=all --output=protein_GDP_fixed.pdb
    
    success "PDB separation completed"
}

################################################################################
# PROTEIN TOPOLOGY
################################################################################
protein_topology() {
    log "=== PROTEIN TOPOLOGY ==="
    
    init_conda
    conda activate md_core 2>/dev/null
    
    log "Running pdb2gmx (AUTO-SELECTING: option 1 = CHARMM36)..."
    echo "1" | gmx pdb2gmx -f protein_GDP_fixed.pdb -o FtsZ_processed.gro -p topol.top -water tip3p -ignh
    
    success "Protein topology created"
}

################################################################################ 
# LIGAND TOPOLOGY
################################################################################
ligand_topology() {
    log "=== LIGAND TOPOLOGY ==="
    
    log "Step 1: Convert to MOL2 in CHIMERA"
    log "  - Open LIG.pdb in CHIMERA"
    log "  - Tools > Structure Editing > AddH"
    log "  - File > Save Mol2 as LIG.mol2"
    read -p "Press Enter when done, or 'q' to skip: " response
    [[ "$response" == "q" ]] && return 1
    
    if [[ ! -f "LIG.mol2" ]]; then
        error "LIG.mol2 not found"
        return 1
    fi
    
    sed -i '' 's/LIG\.pdb/LIG/g' LIG.mol2
    
    log "Step 2: Submit to CGenFF at https://cgenff.com"
    log "  - Upload LIG.mol2, use version 5.0"
    log "  - Download and extract results"
    read -p "Press Enter when done: " response
    
    if [[ ! -f "LIG.str" ]] || [[ ! -f "LIG.cgenff.mol2" ]]; then
        error "CGenFF output files not found"
        return 1
    fi
    
    if [[ ! -f "cgenff_charmm2gmx_py3_nx2.py" ]]; then
        error "cgenff_charmm2gmx_py3_nx2.py not found"
        return 1
    fi
    
    init_conda
    conda activate md_core 2>/dev/null
    
    log "Converting CGenFF to GROMACS format..."
    python3 cgenff_charmm2gmx_py3_nx2.py LIG LIG.cgenff.mol2 LIG.str charmm36-jul2022.ff
    
    success "Ligand topology created"
}

################################################################################
# PREPARE LIGAND COORDINATES
################################################################################
prepare_ligand() {
    log "=== PREPARING LIGAND COORDINATES ==="
    
    init_conda
    conda activate md_core 2>/dev/null
    
    log "Removing virtual sites from lig.itp..."
    if grep -q "vsites" lig.itp; then
        log "  Please manually remove virtual site lines from lig.itp"
        read -p "Press Enter when done: " response
    fi
    
    sed -i '' '/^[[:space:]]*60[[:space:]]/d' lig.itp
    sed -i '' '/^\[ *virtual_sites2 *\]/,$d' lig.itp
    
    log "Converting ligand formats..."
    eval "$(conda shell.zsh hook)" 2>/dev/null || eval "$(conda shell.bash hook)"
    conda activate obabel_env 2>/dev/null
    obabel LIG.cgenff.mol2 -O lig_h.pdb
    
    conda activate md_core 2>/dev/null
    gmx editconf -f lig_h.pdb -o lig_h.gro
    
    log "Creating complex..."
    gmx editconf -f FtsZ_processed.gro -o FtsZ_processed.pdb 
    cat FtsZ_processed.pdb lig_h.pdb > complex_correct.pdb
    log "Please make complex_correct.pdb continuous, only have ATOM lines, remove MODEL, ENDMDL, MASTER, TITLE, REMARK, CRYST1, SEQRES, and CONECT records, etc "
    read -p "Press Enter when done, or 'q' to skip: " response
    [[ "$response" == "q" ]] && return 1
    
    success "complex_correct.pdb ready"

    gmx editconf -f complex_correct.pdb -o complex_correct.gro
    
    success "Ligand coordinates ready"
}

################################################################################
# SETUP SIMULATION
################################################################################
setup_simulation() {
    log "=== SETTING UP SIMULATION ==="
    
    init_conda
    conda activate md_core 2>/dev/null

    log "Creating simulation box..."
    gmx editconf -f complex_correct.gro -o complex_box.gro -bt cubic -d 1.0

    log "Please manually open topol.top, after the charmm36.itp..., type:
         ;Include ligand parameters
         #include \"lig.prm\"
         #include \"lig.itp\" and under [molecules], type
         LIG   1, and save the results"
    read -p "Press Enter when done, or 'q' to skip: " response
    [[ "$response" == "q" ]] && return 1
    success "topol.top ready"
    
    log "Solvating system..."
    gmx solvate -cp complex_box.gro -cs spc216.gro -o solv.gro -p topol.top

    log "Please manually open topol.top, make [molecules], SOL is on new line under LIG, type:"
    read -p "Press Enter when done, or 'q' to skip: " response
    [[ "$response" == "q" ]] && return 1
    success "topol.top ready"

    log "fix the ions.prm, by adding the lines below in the correct category, MAKE SURE INDENT MATCHES
         
        [ bondtypes ]
        ;      i        j  func           b0           kb
          CG2DC1    CG324     1     0.150800    376560.00
          CG2DC1    CG324     1     0.150800    376560.00


        [ angletypes ]
        ;      i        j        k  func       theta0       ktheta          ub0          kub
          CG2R66   CG2R61    CLGR1     5   120.000000   502.080000   0.00000000         0.00
          CG2DC1    CG324    NG3P1     5   110.500000   418.400000   0.24850000     25940.80
          CG2R64    NG311   CG2R61     5   117.000000   418.400000   0.00000000         0.00
          CG2DC1   CG2DC1    CG324     5   123.500000   418.400000   0.00000000         0.00
           CG324   CG2DC1     HGA4     5   116.000000   293.080000   0.00000000         0.00
          CG2DC1    CG324     HGA2     5   109.500000   293.080000   0.00000000         0.00
          CG2R64    NG311   CG2R61     5   117.000000   418.400000   0.00000000         0.00
          CG2DC1   CG2DC1    CG324     5   123.500000   418.400000   0.00000000         0.00
           CG324   CG2DC1     HGA4     5   116.000000   293.080000   0.00000000         0.00
          CG2DC1    CG324     HGA2     5   109.500000   293.080000   0.00000000         0.00



        [ dihedraltypes ]
        ;      i        j        k        l  func         phi0         kphi  mult
           CG2O1   CG2DC1   CG2DC1    CG324     9   180.000000     2.343040     1
           CG2O1   CG2DC1   CG2DC1    CG324     9   180.000000    29.288000     2
          CG2DC1   CG2DC1    CG324    NG3P1     9   180.000000    19.413760     2
          CG2DC1   CG2DC1    CG324    NG3P1     9     0.000000     1.631760     3
            HGA4   CG2DC1    CG324    NG3P1     9   180.000000     6.527040     2
            HGA4   CG2DC1    CG324    NG3P1     9     0.000000     1.506240     3
           NG311   CG2R61   CG2R61    HGR62     9   180.000000    10.041600     2
           CLGR1   CG2R61   CG2R66   CG2R61     9   180.000000    12.552000     2
           CLGR1   CG2R61   CG2R66     FGR1     9   180.000000     8.326160     2
          CG2DC1    CG324    NG3P1    CG324     9     0.000000     0.000000     1
          CG2DC1    CG324    NG3P1    CG324     9   180.000000     8.033280     3
          CG2DC1    CG324    NG3P1     HGP2     9   180.000000     2.175680     3
          CG2R61    NG311   CG2R64   NG2R62     9   180.000000    10.041600     2
          CG2R61    NG311   CG2R64   CG2R61     9   180.000000    10.041600     2
          CG2R64    NG311   CG2R61   CG2R61     9   180.000000    10.041600     2
          CG2R66   CG2R61   CG2R61    HGR62     9   180.000000    10.041600     2
            HGA4   CG2DC1   CG2DC1    CG324     9   180.000000    19.413760     2
            HGA4   CG2DC1   CG2DC1    CG324     9     0.000000     1.631760     3



        [ dihedraltypes ]
        ; 'improper' dihedrals 
        ;      i        j        k        l  func         phi0         kphi
          CG2R64   CG2R61   NG2R62    NG311     2     0.000000   334.720000"

    read -p "Press Enter when done, or 'q' to skip: " response
    [[ "$response" == "q" ]] && return 1
    success "lig.prm ready"

    log "Adding ions..."
    gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -r solv.gro
    log "AUTO-SELECTING: SOL (water molecules for ion replacement)..."
    echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -pname NA -nname CL -neutral
    
    log "Update the topol.top if didn't update, [molecules], remove the # of charges from SOL, and under type
         NA  19"
    read -p "Press Enter when done, or 'q' to skip: " response
    [[ "$response" == "q" ]] && return 1
    success "topol.top ready"

    success "Simulation setup completed"
}

################################################################################
# ENERGY MINIMIZATION
################################################################################
energy_minimization() {
    log "=== ENERGY MINIMIZATION ==="
    
    if [[ ! -f "em.mdp" ]]; then
        error "em.mdp not found"
        return 1
    fi
    
    init_conda
    conda activate md_core 2>/dev/null
    
    log "Running energy minimization..."
    gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
    gmx mdrun -v -deffnm em
    
    log "Extracting potential energy..."
    echo "11 0" | gmx energy -f em.edr -o potential.xvg
    
    success "Energy minimization completed"
}

################################################################################
# RESTRAINING LIGAND
################################################################################
restraint_ligand() {
    log "=== RESTRAINING LIGAND ==="
    
    log "Creating ligand index file..."
    {
        echo "2& ! a H*"
        echo "q"
    } | gmx make_ndx -f lig_h.gro -o index_lig.ndx
    
    log "Generating position restraints..."
    echo "3" | gmx genrestr -f lig_h.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000
    
    log "Adding restraint includes to lig.itp..."
    cat >> lig.itp << 'EOF'

; Ligand position restraints
#ifdef POSRES_LIG
#include "posre_lig.itp"
#endif
EOF
    
    success "Ligand restraints created"
}

################################################################################
# NVT EQUILIBRATION
################################################################################
nvt_equilibration() {
    log "=== NVT EQUILIBRATION ==="
    
    if [[ ! -f "nvt.mdp" ]]; then
        error "nvt.mdp not found"
        return 1
    fi
    
    init_conda
    conda activate md_core 2>/dev/null
    
    log "Running NVT equilibration..."
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
    gmx mdrun -deffnm nvt
    
    log "Extracting temperature data..."
    echo "16 0" | gmx energy -f nvt.edr -o temperature.xvg
    
    success "NVT equilibration completed"
}

################################################################################
# NPT EQUILIBRATION
################################################################################
npt_equilibration() {
    log "=== NPT EQUILIBRATION ==="
    
    if [[ ! -f "npt.mdp" ]]; then
        error "npt.mdp not found"
        return 1
    fi
    
    init_conda
    conda activate md_core 2>/dev/null
    
    log "Running NPT equilibration..."
    gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr -maxwarn 1
    gmx mdrun -deffnm npt
    
    log "Extracting pressure data..."
    echo "17 0" | gmx energy -f npt.edr -o pressure.xvg
    
    log "Extracting density data..."
    echo "23 0" | gmx energy -f npt.edr -o density.xvg
    
    success "NPT equilibration completed"
}

################################################################################
# MD PRODUCTION
################################################################################
md_production() {
    log "=== MD PRODUCTION ==="
    
    if [[ ! -f "md.mdp" ]]; then
        error "md.mdp not found"
        return 1
    fi
    
    init_conda
    conda activate md_core 2>/dev/null
    
    log "Running MD production (this may take a while)..."
    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1
    gmx mdrun -deffnm md
    
    success "MD production completed"
}

################################################################################
# FIX TRAJECTORY
################################################################################
fix_trajectory() {
    log "=== FIXING TRAJECTORY ==="
    
    init_conda
    conda activate md_core 2>/dev/null
    
    log "Centering and removing PBC (AUTO-SELECTING: Protein, System)..."
    echo "1 0" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol -ur compact
    
    log "Extracting final frame (AUTO-SELECTING: System)..."
    echo "0" | gmx trjconv -s md.tpr -f md_center.xtc -o final_frame.pdb -dump 1
    
    log "Fitting trajectory (AUTO-SELECTING: Backbone, System)..."
    echo "4 0" | gmx trjconv -s md.tpr -f md_center.xtc -o md_fit.xtc -fit rot+trans
    
    success "Trajectory files fixed"
}

################################################################################
# POST-PROCESSING
################################################################################
postprocessing() {
    log "=== POST-PROCESSING ANALYSIS ==="
    
    init_conda
    conda activate md_core 2>/dev/null
    
    log "Step 1: Creating index file for ligand heavy atoms..."
    {
        echo "r LIG & ! a H*"
        echo "name 21 LIG_heavy"
        echo "q"
    } | gmx make_ndx -f md.tpr -o index.ndx
    
    log "Step 2: Protein backbone RMSD (AUTO-SELECTING: Protein, Protein)..."
    echo "1 1" | gmx rms -s md.tpr -f md_fit.xtc -n index.ndx -o rmsd_FtsZ.xvg
    
    log "Step 3: Ligand RMSD (AUTO-SELECTING: Protein, LIG_heavy)..."
    echo "1 21" | gmx rms -s md.tpr -f md_fit.xtc -n index.ndx -o rmsd_LIG.xvg
    
    log "Step 4: Hydrogen bonds..."
    gmx hbond -s md.tpr -f md_fit.xtc -n index.ndx -r Protein -t LIG -num hbond_protein_lig.xvg
    
    log "Step 5: Radius of gyration (AUTO-SELECTING: Protein)..."
    echo "1" | gmx gyrate -s md.tpr -f md_fit.xtc -o gyrate.xvg
    
    log "Step 6: Contacts/Packing (AUTO-SELECTING: LIG_heavy)..."
    echo "21" | gmx mindist -s md.tpr -f md_fit.xtc -n index.ndx
    
    log "Step 7: SASA (AUTO-SELECTING: LIG_heavy)..."
    echo "21" | gmx sasa -s md.tpr -f md_fit.xtc -n index.ndx
    mv area.xvg sasa.xvg
    
    log "Step 8: Protein flexibility RMSF (AUTO-SELECTING: Protein)..."
    echo "1" | gmx rmsf -s md.tpr -f md_fit.xtc -n index.ndx -res
    
    log "Step 9: Cluster analysis (AUTO-SELECTING: LIG_heavy)..."
    echo "21" | gmx cluster -s md.tpr -f md_fit.xtc -n index.ndx
    
    success "Post-processing analysis completed"
}

################################################################################
# PLOTTING
################################################################################
plot_results() {
    log "=== GENERATING PLOTS ==="
    
    if [[ ! -f "plot.py" ]]; then
        error "plot.py not found"
        return 1
    fi
    
    init_conda
    eval "$(conda shell.zsh hook)" 2>/dev/null || eval "$(conda shell.bash hook)"
    conda activate md_analysis 2>/dev/null
    
    log "Running plotting script..."
    python plot.py
    
    success "Plots generated in Figures/ directory"
}

################################################################################
# MAIN MENU
################################################################################
main_menu() {
    echo ""
    echo "=========================================="
    echo "MOLECULAR DYNAMICS WORKFLOW - Mac Version"
    echo "        (AUTO-SELECTING GROUP OPTIONS)"
    echo "=========================================="
    echo ""
    echo "Select workflow:"
    echo "1. Check prerequisites only"
    echo "2. PDB preparation and separation"
    echo "3. Create topologies (protein + ligand)"
    echo "4. Prepare ligand coordinates"
    echo "5. Setup simulation (box, solvation, ions)"
    echo "6. Energy minimization"
    echo "7. Equilibration (NVT + NPT)"
    echo "8. MD production run"
    echo "9. Analysis (RMSD, H-bonds, etc)"
    echo "10. Full workflow (all steps)"
    echo "11. Exit"
    echo ""
    read -p "Enter choice (1-11): " choice
    
    case $choice in
        1) check_prerequisites ;;
        2) 
            check_prerequisites
            prepare_pdb
            separate_pdb
            ;;
        3)
            check_prerequisites
            protein_topology
            ligand_topology
            ;;
        4)
            check_prerequisites
            prepare_ligand
            ;;
        5)
            check_prerequisites
            setup_simulation
            ;;
        6)
            check_prerequisites
            energy_minimization
            restraint_ligand
            ;;
        7)
            check_prerequisites
            nvt_equilibration
            npt_equilibration
            ;;
        8)
            check_prerequisites
            md_production
            ;;
        9)
            check_prerequisites
            fix_trajectory
            postprocessing
            ;;
        10)
            check_prerequisites
            prepare_pdb
            separate_pdb
            protein_topology
            ligand_topology
            prepare_ligand
            setup_simulation
            energy_minimization
            restraint_ligand
            nvt_equilibration
            npt_equilibration
            md_production
            fix_trajectory
            postprocessing
            plot_results
            ;;
        11) exit 0 ;;
        *)
            error "Invalid choice"
            main_menu
            ;;
    esac
    
    success "✅ Completed!"
    echo ""
    main_menu
}

# Start script
main_menu
