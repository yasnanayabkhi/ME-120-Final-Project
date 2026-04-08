#!/bin/bash

#---------------------------------------------------------------
# Script for ionization of Dacomitinib-FtsZ complex
# 3rd step of MD: ionization
# Files required: ions.mdp, FtsZ-Dacomitinib_solv.gro, topol.top
#---------------------------------------------------------------

# To create atomic description of FtsZ-Dacomitinib system
# Generates 1 file: ions.tpr 
gmx grompp ions.mdp -c FtsZ-Dacomitinib_solv.gro -p topol.top -o ions.tpr 

# Pass file to genion & adds ions (sodium, chloride) to system
# Generates 1 file: FtsZ-Dacomitinib_solv_ions.gro 
# echo 13 for group 13 "SOL" to replace small amount of sol with ions instead to balance charge
echo 13 | gmx genion -s ions.tpr -o FtsZ-Dacomitinib_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# tell user ionization is done 
echo "ionization is completed.\n"
echo "The file generated and edited is ions.tpr & FtsZ-Dacomitinib_solv_ions.gro"
