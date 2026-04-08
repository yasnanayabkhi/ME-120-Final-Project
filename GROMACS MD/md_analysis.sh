#!/bin/bash

#----------------------------------------------------
# Script for MD analysis of Dacomitinib-FtsZ system
# 7th step of MD: MD analysis 
# Files required: md.tpr, md_noPBC.xtc, em.tpr
#----------------------------------------------------

# To calculate Root-Mean-Square Deviation (RMSD) for backbone & crystal structure
# Generates 1 file: rmsd.xvg, rmsd_xtal.xvg
# echo "4\n4" , group 4 is Backbone for least-square fit & calculation
echo "4\n4" | gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns
gmx rms -s em.tpr -f md_noPBC.xtc -o rmsd_xtal.xvg -tu ns

# To calculate Root-Mean-Square Fluctuations (RMSF) according to Alpha Carbon
# Generates 1 file: rmsf.xvg
# echo 6 , group 6 is alpha carbon 
echo 6 | gmx rmsf -s md.tpr -f md_noPBC.xtc -o rmsf.xvg, -res

# To calculate Radius of Gyration (analyzing compactness)
# Generates 1 file: gyrate.xvg
gmx gyrate -s md.tpr -f md_noPBC.xtc -o gyrate.xvg -sel Protein -tu ns

# To determine the secondary structure count during simulation & calculate DSSP
# Generates 2 files: dssp.dat, dssp_num.xvg
gmx dssp -s md.tpr -f md_noPBC.xtc -tu ns -o dssp.dat -num dssp_num.xvg

# To determine hydrogen bonding through the simulation
# Generates 3 files: hbnum_mainchain.xvg, hbnum_sidechain.xvg, hbnum_prot_wat.xvg
# echo "7\n7", group 7 (MainChain+H) for both
# echo "8\n8", group 8 (SideChain) for both
# echo "1\n12", group 1 (Protein) and group 12 (Water)
echo "7\n7" | gmx hbond -s md.tpr -f md_noPBC.xtc -tu ns -num hbnum_mainchain.xvg
echo "8\n8" | gmx hbond -s md.tpr -f md_noPBC.xtc -tu ns -num hbnum_sidechain.xvg
echo "1\n12" | gmx hbond -s md.tpr -f md_noPBC.xtc -tu ns -num hbnum_prot_wat.xvg

# tell user MD analysis is done 
echo "MD analysis is completed.\n"
echo "The files generated are rmsd.xvg, rmsd_xtal.xvg, rmsf.xvg, gyrate.xvg, dssp.dat, dssp_num.xvg, hbnum_mainchain.xvg, hbnum_sidechain.xvg, hbnum_prot_wat.xvg"
