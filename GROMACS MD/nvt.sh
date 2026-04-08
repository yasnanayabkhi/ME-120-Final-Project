#!/bin/bash

#-------------------------------------------
# Script for NVT of Dacomitinib-FtsZ system
# 5th step of MD: NVT
# Files required: nvt.mdp, em.gro
#-------------------------------------------

# To create nvt.tpr
# Generates 1 file: nvt.tpr
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

# To run the NVT equilibration (310K)
# Generates 5 files: nvt.log, nvt.edr, nvt.trr, nvt.gro, nvt.cpt
gmx mdrun -deffnm nvt

# To analyze the temperature progress during NVT steps
# Generates 1 file: temperature.xvg
# echo 16 0 to choose group 16 (temperature) then 0 to exit
# (can be plotted using md_figures.ipynb if needed)
gmx energy -f nvt.edr -o temperature.xvg

# tell user NVT equilibration is done 
echo "NVT equilibration is completed.\n"
echo "The file generated are nvt.tpr, nvt.log, nvt.edr, nvt.trr, nvt.gro, nvt.cpt, temperature.xvg"