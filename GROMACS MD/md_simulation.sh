#!/bin/bash

#----------------------------------------------------
# Script for MD simulation of Dacomitinib-FtsZ system
# 7th step of MD: MD simulation
# Files required: md.mdp, npt.gro, npt.cpt, topol.top
#----------------------------------------------------

# To create md.tpr
# Generates 1 file: md.tpr
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

# To run the MD simulation (10 ns at 310K)
# Generates 5 files: md.xtc, md.cpt, md.edr, md.gro, md.log
gmx mdrun -deffnm md

# Correct for Periodicity Effects
# Generates 1 file: md_noPBC.xtc
gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center

# tell user MD simulation is done 
echo "MD simulation is completed.\n"
echo "The file generated are md.tpr, md.xtc, md_noPBC.xtc, md.cpt, md.edr, md.gro, md.log"
