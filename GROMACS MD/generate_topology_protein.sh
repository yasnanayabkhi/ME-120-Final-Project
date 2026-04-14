#!/bin/bash

#--------------------------------------------------------
# Script for generating topology of FtsZ
# 1st step of MD: generate FtsZ topology
# Files required: FtsZ-Dacomitinib128.pdb
#--------------------------------------------------------

# To save seperate coordinates to isolate dacomitinib from complex 
# Generates 1 file: dacomitinib.pdb
grep LIG FtsZ-Dacomitinib128.pdb > dacomitinib.pdb

# To remove LIG (Dacomitinib) lines and isolate new file for just FtsZ
# Generates 1 file: FtsZ.pdb
grep -v "LIG" FtsZ-Dacomitinib128.pdb > FtsZ.pdb

# To fix ASN 25 residue (before doing pdb2gmx) because missing atom CG
pdbfixer FtsZ.pdb --add-atoms=all --output=FtsZ_fixed.pdb

# Generates the topology for FtsZ
# Generates 1 file: FtsZ_processed.gro and topol.top (FtsZ)
# echo 1 for charmm36 force field
# -ignh: ignores hydrogen then rebuilds using charmm36 rules
echo "1" | gmx pdb2gmx -f FtsZ_fixed.pdb -o FtsZ_processed.gro -water tip3p -ignh

# tell user generating protein topology is done 
echo "FtsZ Topology generation is completed.\n"
echo "The files generated are dacomitinib.pdb, FtsZ.pdb, FtsZ_processed.gro, and topol.top"
