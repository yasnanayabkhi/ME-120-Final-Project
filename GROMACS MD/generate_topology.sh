#!/bin/bash

#--------------------------------------------------------
# Script for generating topology of Dacomitinib-FtsZ
# 1st step of MD: generate topology
# Files required: FtsZ-Dacomitinib.pdb
#--------------------------------------------------------

# Clean the PDB (removal of water)
# Generates 3 files: topology file, position restraint file, post-procssed structure file
grep -v HOH FtsZ-Dacomitinib128.pdb > FtsZ-Dacomitinib_clean.pdb

# To fix ASN 25 residue (before doing pdb2gmx) 
pdbfixer FtsZ-Dacomitinib_clean.pdb --output=FtsZ-Dacomitinib_fixed.pdb --add-atoms=all --add-residues

# Processed with tip3p water
# Generates 1 file: processed
# echo 1 to choose 1st choice for force field (charmm36)
echo 1 | gmx pdb2gmx -f FtsZ-Dacomitinib_fixed.pdb -o FtsZ-Dacomitinib_processed.gro -water tip3p

# tell user generating topology is done 
echo "Topology generation is completed.\n"
echo "The files generated are FtsZ-Dacomitinib_clean.pdb and FtsZ-Dacomitinib_processed.gro"
