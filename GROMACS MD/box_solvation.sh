#!/bin/bash

#----------------------------------------------
# Script to create box and solvate
# 2nd step of MD: Define box & Solvate
# Files required: FtsZ-Dacomitinib_procssed.gro
#-----------------------------------------------

# Create the box of 1.2nm buffer & cubic shape
# Generates 1 file: FtsZ-Dacomitinib_newbox.gro
gmx editconf -f FtsZ-Dacomitinib_processed.gro -o FtsZ-Dacomitinib_newbox.gro -c -d 1.2 -bt cubic

# Solvate the box 
# Generates 1 file: FtsZ-Dacomitinib_solv.gro
gmx solvate -cp FtsZ-Dacomitinib_newbox.gro -cs spc216.gro -o FtsZ-Dacomitinib_solv.gro -p topol.top

