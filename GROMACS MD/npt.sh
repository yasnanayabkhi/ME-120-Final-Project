#!/bin/bash

#-------------------------------------------
# Script for NPT of Dacomitinib-FtsZ system
# 6th step of MD: NPT
# Files required: npt.mdp, nvt.gro
#-------------------------------------------

# To create npt.tpr
# Generates 1 file: npt.tpr
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

# To run the NPT equilibration
# Generates 5 files: npt.log, npt.edr, npt.trr, npt.gro, npt.cpt
gmx mdrun -deffnm nvt

# To analyze the pressure progress during NPT steps
# Generates 1 file: pressure.xvg
# echo 17 0 to choose group 17 (pressure) then 0 to exit
# (can be plotted using md_figures.ipynb if needed)
gmx energy -f npt.edr -o pressure.xvg

# To analyze the density progress during NPT steps
# Generates 1 file: density.xvg
# echo 23 0 to choose group 23 (density) then 0 to exit
# (can be plotted using md_figures.ipynb if needed)
gmx energy -f npt.edr -o density.xvg

# tell user NPT equilibration is done 
echo "NPT equilibration is completed.\n"
echo "The file generated are npt.tpr, npt.log, npt.edr, npt.trr, npt.gro, npt.cpt, temperature.xvg, density.xvg"