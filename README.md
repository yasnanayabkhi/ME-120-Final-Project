# Computational Study of Dacomitinib Inhibition of Staphylococcus Aureus FtsZ using Molecular Docking and Molecular Dynamics


## Project Summary

Problem: 
This project will investigate the inhibition potential of Dacomitinib by understanding the dynamic stability of 
Dacomitinib-FtsZ complex. FtsZ is an essential protein for cell division in S. Aureus, including methicillin-resistant strains (MRSA). Inhibiting FtsZ can essentially stop bacterial proliferation which becomes a potential therapeutic target for developing a new antibiotic. While Dacomitinib is an approved cancer drug, virtual screening & experimental studies have suggested that Dacomitinib does inhibit FtsZ. The dynamic behavior of this protein-ligand complex has not been characterized before. 

Computation Methods: Molecular docking (Webina - Autodock Vina) and Molecular Dynamics (GROMACS) to model & analyze FtsZ-Dacomitinib complex interactions over 1 nanosecond simulation. 

Main Results: The main outputs are the static binding affinity (e.g. stdout128.txt), docked complex (FtsZ-Dacomitinib128.pdb), simulation trajectories (md.xtc, md_center.xtc, and md_fit.xtc), quantitative analyses (e.g. RMSF, RMSD, hydrogen bonding, interaction energy, etc.), and the generated plots (e.g. Dacomitinib_FtsZ_mindist.png, etc.) to assess the binding stability of the ligand to the protein. This could potentially provide insight into the potential of Dacomitinib as an inhibitor to S. Aureus FtsZ.


## Package Contents

    - README.md: documentation describing workflow, setup, and reproducibility instructions

    - Cleaning/: directory of workflow for setting up for docking and includes cleaned PDB and SDF as well. 

    - Docking/: directory of workflow for obtaining docked complex to use for MD

    - GROMACS MD/: directory of workflow for obtaining 1 ns MD simulation and generated plots

    - GROMACS MD/yml_files/: directory that includes all .yml files for the 3 conda enviornments

    - md_core.yml: conda environment file for MD simulation workflow (includes required tools with GROMACS)

    - obabel_env.yml: conda environment file for ligand conversion (mol2 to pdb)

    - md_analysis.yml: conda environment file for MD figure generation (paired with plot.py w/ needed dependencies)

    - GROMACS MD/mdp_files/: directory that includes all .mdp files for the steps of MD (e.g. md.mdp, em.mdp, etc.)

    - md_workflow.sh: automated script to run full molecular dynamics pipeline (setup --> simulation --> processing (plots))

    - plot.py: python script for post-processing of .xvg files and generating figures such as RMSD, RMSF, etc.

    - cgenff_charmm2gmx_py3_nx2.py: to convert LIG cgenff files into gromacs format to be usable for MD (generated lip.itp and lig.prm)

    - charmm36-jul2022.ff/: force field directory used for MD simulation

    - Results/: directory that contains all 3 runs with the results from docking and MD including the generated plots.

    - Results/result_""/MD/figures: directory containing generated plots via plot.py 

    - *.xtc: compressed trajectory files from the 1 ns MD simulation

    - *.pdb, *.gro: structure files used at different stages of the workflow

    - index.ndx: index ile that contains defined atom groups for analysis (includes LIG_heavy group)

    - topol.top: main topology file describing force field parameters & the system

    - lig.top: ligand topology generated using cgenff_charmm2gmx_py3_nx2.py 

    - lig.prm: ligand parameter file generated from CGenFF and cgenff_charmm2gmx_py3_nx2.py

    - posre_lig.itp: ligand position restraints file to keeps the ligand from drifting


## Environment Setup

Make sure you have miniconda 3 downloaded:
    - If you don't have it installed, follow procedure on this website: 
    https://www.anaconda.com/docs/getting-started/miniconda/install/overview 

To obtain charmm36 force field:
    - http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs 
        --> download charmm36-jul2022.ff.tgz 
        --> unzip it manually or through terminal, type: tar -zxvf charmm36-jul2022.ff.tgz
        --> save it as charmm36-jul2022.ff into your folder

Download a text editor such as VSCode:
    - https://code.visualstudio.com/download 

Needed 3 conda environments for this MD procedure:
 - Download obabel_env.yml, in terminal, type: 
            conda env create -f obabel_env.yml
- Download md_core.yml, in terminal, type: 
            conda env create -f md_core.yml
- Download md_analysis.yml, in terminal, type: 
            conda env create -f md_analysis.yml 

Download chimera for editing and conversion of pdbs: 
    - https://www.cgl.ucsf.edu/chimera/download.html 

Make an CGenFF account on the website below:
    - https://cgenff.com 


## How To Run

1. In the Cleaning folder, open procedure_for_cleaning.txt

    - Follow this procedure to obtained cleaned files for Docking
    - specifications are given within the procedure

2. In the Docking folder, open procedure_for_docking.txt 

    - Follow this procedure to obtain the files for MD. 
    - specifications are given within the procedure

3. In the GROMACS MD, open procedure_for_md.txt

    - Follow this procedure to obtain the figures from the 1 ns MD simulation
    - Either do automatically using md_workflow.sh or manually using the   procedure below the automatic part in the procedure_for_md.txt
    - specifications are given within the procedure


## Expected Outputs

The main result files are: 
- stdout128.txt (AutoDock Vina Results)
- FtsZ-Dacomitinib128.pdb (Docked Complex using Webina ouput)
- md.gro 
- md.xtc
- md.tpr
- md_center.xtc
- md_fit.xtc 

The main figures are:
- Dacomitinib_FtsZ_gyrate.png
- Dacomitinib_FtsZ_hbond_protein_lig.png 
- Dacomitinib_FtsZ_interaction_energy.png
- Dacomitinib_FtsZ_mindist.png
- Dacomitinib_FtsZ_rmsd_dist.png
- Dacomitinib_FtsZ_rmsd.png
- Dacomitinib_FtsZ_rmsf.png
- Dacomitinib_FtsZ_sasa.png


## Runtime Notes And Limitations

Cleaning: ~5-15 minutes
    - Protein Data Bank: https://www.rcsb.org/ligand/1C9
                         https://www.rcsb.org/structure/4DXD 
    - PyMOL: https://www.pymol.org 


Docking: ~10 minutes to 2 hours (depending on system performance and server load)
    - Webina (Autodock Vina): https://durrantlab.pitt.edu/webina/ 

GROMACS MD: ~5 hours to 2 days (depending on system performance)
    - required conda environments: obabel_env 
                                   md_core
                                   md_analysis 

    - required scripts: md_workflow.sh
                        plot.py 
                        cgenff_charmm2gmx_py3_nx2.py 

    - CGenFF: https://cgenff.com 

Limitations: Workflow requires multiple special dependencies, excellent file organization, and external websites. Multiple steps require manual edits to the files, which can introduce human error. Full reproducibility requires all the scripts, environments, and force field files to be installed correctly to work. High computational cost for MD runs & large storage requirements for trajectory files (.trr ), can reach to high MB size. Requirement of user intervention in editing lig.prm, topol.top, and formatting can introduce potential variability.


## Omitted Files, If Any

- nvt.trr and npt.trr were omitted.

- nvt.trr and npt.trr were omitted due to the large size of the  files, ~37 MB and ~890 MB respectively. These are the full, high-precision simulation trajectory files (alot of data).

- The omission doesn't affect the reproducibility by how all the simulations can be generated fully with use of files such as .tpr, .mdp, .top, and many more files. The results of the MD simulation are still avaliable in the form of .xtc such as md.xtc, md_fit.xtc, etc.


## AI Usage Appendix

# Use 1
- Tool - ChatGPT-4o via chat.openai.com, April 2026
- Purpose - Used to fix an outdated cgenff_charmm2gmx_py3_nx2.py for  creating lig.itp, lig.top, lig.pdb and lig.prm from a CGenFF updated lig.mol2 file
- Prompts - “I am using this python script in order to convert a lig.mol2 file into a lig.itp file through CGenFF. Fix this outdated python script to allow for varying versions of Networkx and any other conflicting version errors”
- Output - “Here is the updated script with fixes for NetworkX version compatibility and other common version-related issues:” Then gave me the fixed python script.
- Modifications - Did not modify output. I used the resulting cgenff_charmm2gmx_py3_nx2.py instead of the older version.
- Verification - Tested code and visualized produced pdb, whilst also looking at its contents of the parameter files to see if it generated everything properly. Some modifications to the lig.itp and lig.top were required after use of this script.

# Use 2
- Tool - Claude Haiku 4.5 via claude.ai, April 2026
- Purpose - Used to generate plots based on all the .xvg files made beforehand for data analysis visualization 
- Prompts - “I want you to generate me a script to generate plots based on files given here. I made an index.ndx that gave group 21 LIG_heavy, protein stability was given by rmsd_FtsZ.xvg, ligand binding stability was given by rmsd_LIG.xvg, protein-ligand interaction energy was given by interaction_energy.xvg, hydrogen bonds given by hbond_protein_lig.xvg, radius of gyration was given by gyrate.xvg, contacts/packing was given by gmx mindist -s md.tpr -f md_fit.xtc -n index.ndx, SASA (burial) was given by sasa.xvg, RMSF was given by rmsf.xvg, cluster analysis was given by rmsd-clus.xpm and rmsd-dis.xvg, also did center of mass, I have attached all the files for reference”
- Output - “Here is the script for generating all the plots from the protein-ligand MD simulation, it is called plot.py.
- Modifications - Did not modify output. It had given proper graphs usable for the project. 
- Verification - Tested code by running the plot.py and made sure the generated plots & means made sense for the data when looking into the .xvg files. Double checked the numbers had made sense. 

# Use 3
- Tool - Claude Haiku 4.5 via claude.ai, April 2026
- Purpose - Used to generate a md_workflow.sh to run the entire MD to plot generation into one script for reproducibility specifically for this project.  
- Prompts - “I have attached my procedure_for_md.txt, generate me a script based on this information, make sure it chooses its input automatically based on what I wrote in the procedure”
- Output - “I’ll create a comprehensive .sh script that converts this procedure into an automated workflow. I have built-in pauses for user interaction. Here is the script, md_workflow.sh.”
- Modifications - I had to make some modifications. I had to add more pre-built pauses (thru log) that were user interactions because of the user’s need to manually change some files such as topol.top, lig.prm, etc. For example:
    log "Please manually open topol.top, after the charmm36.itp..., type:
        ;Include ligand parameters
        #include \"lig.prm\"
        #include \"lig.itp\" and under [molecules], type
        LIG   1, and save the results"
   read -p "Press Enter when done, or 'q' to skip: " response
   [[ "$response" == "q" ]] && return 1
   success "topol.top ready"

- Verification - Tested code by running the md_workflow.sh with the input files and seeing it was running properly, if not, I would edit the script to make it work properly.

