#!/bin/bash

#-----------------------------------------------------------
# Script for Energy Minimization of Dacomitinib-FtsZ system
# 4rd step of MD: Energy Minimization
# Files required: minim.mdp, FtsZ-Dacomitinib_solv_ions.gro
#-----------------------------------------------------------

# To create em.tpr 
# Generates 1 file: em.tpr
gmx grompp -f minim.mdp -c FtsZ-Dacomitinib_solv_ions.gro -p topol.top -o em.tpr

# To run the energy minimizaiton 
# Generates 4 files: em.log, em.edr, em.trr, em.gro
gmx mdrun -v -deffnm em

# To analyze the potential energy distribution thru the steps
# Generates 1 file: potential.xvg  
# echo 11 0 to choose group 11 (potential) then 0 terminates input
# (can be plotted using md_figures.ipynb if needed)
echo "11 0" | gmx energy -f em.edr -o potential.xvg 

# tell user energy minimization is done 
echo "energy minimization is completed.\n"
echo "The file generated are em.tpr, em.log, em.edr, em.trr, em.gro, potential.xvg"