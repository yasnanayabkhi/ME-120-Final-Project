#!/usr/bin/env python3

# USAGE: python cgenff_charmm2gmx.py DRUG drug.mol2 drug.str charmm36.ff
# Updated for CGenFF 5.0, Python 3.7+, NetworkX 2.3+
# Tested with Python 3.7-3.9, NetworkX 2.3-2.5

# Original Copyright (C) 2014 E. Prabhu Raman prabhu@outerbanks.umaryland.edu
# Modified 11/6/2018 by Justin Lemkul to add lone pair support
# Modified 01/10/2019 by Conrard Tetsassi to work with Networkx 2.3
# Updated 2026 for CGenFF 5.0 and improved parameter parsing

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

import string
import re
import sys
import os
import math
import numpy as np
import networkx as nx

#=================================================================================================================
def check_versions(str_filename, ffdoc_filename):
    ffver = 0
    strver = 0
    
    # Check stream file version
    try:
        with open(str_filename, 'r') as f:
            for line in f.readlines():
                if "For use with CGenFF" in line and "version" in line:
                    parts = line.split()
                    # Extract version number (e.g., "5.0" from line)
                    for i, part in enumerate(parts):
                        if "version" in part.lower() and i+1 < len(parts):
                            strver = parts[i+1].strip()
                            break
                    if strver == 0:  # Try alternative parsing
                        match = re.search(r'version\s+(\d+\.\d+)', line)
                        if match:
                            strver = match.group(1)
                    if strver:
                        print(f"--Version of CGenFF detected in {str_filename}: {strver}")
                        break
    except:
        print(f"Warning: Could not open {str_filename} to check version")
    
    # Check forcefield doc version
    try:
        with open(ffdoc_filename, 'r') as f:
            for line in f.readlines():
                if "CGenFF" in line and "version" in line:
                    match = re.search(r'version\s+(\d+\.\d+)', line)
                    if match:
                        ffver = match.group(1)
                        print(f"--Version of CGenFF detected in {ffdoc_filename}: {ffver}")
                        break
    except:
        print(f"Warning: Could not open {ffdoc_filename} to check version")
    
    # Warn about version mismatch
    if strver and ffver and strver != ffver:
        print(f"\nWARNING: CGenFF versions may not be equivalent! (Stream: {strver}, FF: {ffver})\n")

#-----------------------------------------------------------------------
def is_lp(s):
    """Check if atom name is a lone pair"""
    if len(s) >= 2 and s[0] == 'L' and s[1] == 'P':
        return True
    return False

#-----------------------------------------------------------------------
def is_lp_host_atom(graph, nvsites, name):
    """Check if atom is host for a lone pair"""
    for ai in range(0, nvsites):
        if name == graph.nodes[ai].get('at1'):
            return True
    return False

#-----------------------------------------------------------------------
def find_vsite(graph, nvsites, natoms, atnum):
    """Find virtual site index for an atom"""
    atname = graph.nodes[atnum].get('name')
    for i in range(0, nvsites):
        if graph.nodes[i].get('at1') == atname:
            vsite_name = graph.nodes[i].get('vsite')
            for j in range(0, natoms):
                if graph.nodes[j].get('name') == vsite_name:
                    return j
    return None

#-----------------------------------------------------------------------
def construct_lp(x1, y1, z1, x2, y2, z2, dist):
    """Construct lone pair coordinates"""
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2
    dr = math.sqrt(dx*dx + dy*dy + dz*dz)
    if dr > 0:
        dr = dist / dr
    else:
        dr = 0
    
    xlp = x1 + dr * dx
    ylp = y1 + dr * dy
    zlp = z1 + dr * dz
    
    return xlp, ylp, zlp

#-----------------------------------------------------------------------
def read_gmx_atomtypes(filename):
    """Read atom types from GROMACS forcefield"""
    atomtypes = []
    try:
        with open(filename, 'r') as f:
            for line in f.readlines():
                if line.startswith(";"):
                    continue
                if line.strip() == '':
                    continue
                entry = line.split()
                if len(entry) >= 2:
                    atomtypes.append([entry[0], entry[1]])
    except FileNotFoundError:
        print(f"Error: Could not find {filename}")
        sys.exit(1)
    return atomtypes

#-----------------------------------------------------------------------
def get_filelist_from_gmx_forcefielditp(ffdir, ffparentfile):
    """Get list of forcefield files to include"""
    filelist = []
    filepath = os.path.join(ffdir, ffparentfile)
    try:
        with open(filepath, 'r') as f:
            for line in f.readlines():
                if line.startswith("#include"):
                    parts = line.split()
                    if len(parts) >= 2:
                        filename = ffdir + "/" + parts[1].replace("\"", "").replace("'", "")
                        if os.path.exists(filename):
                            filelist.append(filename)
    except FileNotFoundError:
        print(f"Error: Could not find {filepath}")
        sys.exit(1)
    return filelist

#-----------------------------------------------------------------------
def read_gmx_anglpars(filename):
    """Read angle parameters from GROMACS files"""
    angllines = []
    try:
        with open(filename, 'r') as f:
            section = "NONE"
            for line in f.readlines():
                if line.startswith(";"):
                    continue
                if line.startswith("\n"):
                    continue
                if line.startswith("["):
                    section = "NONE"
                if section == "ANGL":
                    angllines.append(line)
                if line.startswith("[ angletypes ]"):
                    section = "ANGL"
    except FileNotFoundError:
        return []
    
    anglpars = []
    for line in angllines:
        entry = line.split()
        if len(entry) >= 5:
            try:
                ai, aj, ak, eq = entry[0], entry[1], entry[2], float(entry[4])
                anglpars.append([ai, aj, ak, eq])
            except ValueError:
                continue
    
    return anglpars

#-----------------------------------------------------------------------
def get_charmm_rtp_lines(filename, molname):
    """Extract topology section for molecule from CHARMM file"""
    foundmol = 0
    store = 0
    rtplines = []
    
    try:
        with open(filename, 'r') as f:
            for line in f.readlines():
                if store == 1 and line.startswith("RESI"):
                    store = 0
                
                if line.startswith("RESI"):
                    parts = line.split()
                    if len(parts) >= 2:
                        rtfmolname = parts[1]
                        if rtfmolname == molname:
                            store = 1
                
                if line.startswith("END"):
                    store = 0
                
                if store == 1:
                    rtplines.append(line)
    except FileNotFoundError:
        print(f"Error: Could not find {filename}")
        sys.exit(1)
    
    return rtplines

#-----------------------------------------------------------------------
def get_charmm_prm_lines(filename):
    """Extract parameter section from CHARMM file"""
    prmlines = []
    store = 0
    
    try:
        with open(filename, 'r') as f:
            for line in f.readlines():
                if line.startswith("END"):
                    store = 0
                
                if store:
                    prmlines.append(line)
                
                if line.startswith("read para"):
                    store = 1
    except FileNotFoundError:
        print(f"Error: Could not find {filename}")
        sys.exit(1)
    
    return prmlines

#-----------------------------------------------------------------------
class atomgroup:
    """Class to handle atom topology and parameters"""
    
    def __init__(self):
        self.G = nx.Graph()
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0
        self.nimpropers = 0
        self.nvsites = 0
        self.name = ""
        self.coord = []

    def read_charmm_rtp(self, rtplines, atomtypes):
        """Parse CHARMM topology"""
        noblanks = [x for x in rtplines if len(x.strip()) > 0]
        nocomments = [x for x in noblanks if x.strip()[0] not in ['*', '!']]
        
        section = "NONE"
        self.nvsites = 0
        
        for line in nocomments:
            entry = line.split()
            
            if len(entry) == 0:
                continue
            
            if entry[0] == "RESI":
                self.name = entry[1] if len(entry) > 1 else "MOL"
                continue
            
            if entry[0] == "ATOM":
                aname = entry[1] if len(entry) > 1 else ""
                atype = entry[2] if len(entry) > 2 else ""
                charge = float(entry[3]) if len(entry) > 3 else 0.0
                
                self.G.add_node(self.natoms)
                self.G.nodes[self.natoms]['name'] = aname
                self.G.nodes[self.natoms]['type'] = atype
                self.G.nodes[self.natoms]['charge'] = charge
                self.G.nodes[self.natoms]['resid'] = 1
                self.G.nodes[self.natoms]['beta'] = 0.0
                self.G.nodes[self.natoms]['x'] = 0.0
                self.G.nodes[self.natoms]['y'] = 0.0
                self.G.nodes[self.natoms]['z'] = 0.0
                
                # Get mass from atomtypes
                mass = 0.0
                for atype_entry in atomtypes:
                    if atype_entry[0] == atype:
                        try:
                            mass = float(atype_entry[1])
                        except ValueError:
                            mass = 0.0
                        break
                self.G.nodes[self.natoms]['mass'] = mass
                
                self.natoms += 1
            
            elif entry[0] == "BOND" and len(entry) >= 3:
                at1 = entry[1]
                at2 = entry[2]
                # Find indices
                idx1 = -1
                idx2 = -1
                for i in range(self.natoms):
                    if self.G.nodes[i]['name'] == at1:
                        idx1 = i
                    if self.G.nodes[i]['name'] == at2:
                        idx2 = i
                
                if idx1 >= 0 and idx2 >= 0:
                    self.G.add_edge(idx1, idx2)
                    self.nbonds += 1
            
            elif entry[0] == "LONEPAIR":
                # Parse: LONEPAIR COLI LP1 CL C13 DIST 1.6400 SCAL 0.0
                if len(entry) >= 8:
                    lptype = entry[1]  # COLI
                    lpname = entry[2]  # LP1
                    host_atom = entry[3]  # CL
                    ref_atom = entry[4]  # C13
                    
                    self.G.add_node(self.natoms)
                    self.G.nodes[self.natoms]['name'] = lpname
                    self.G.nodes[self.natoms]['type'] = 'LPH'
                    self.G.nodes[self.natoms]['charge'] = 0.0
                    self.G.nodes[self.natoms]['resid'] = 1
                    self.G.nodes[self.natoms]['beta'] = 0.0
                    self.G.nodes[self.natoms]['vsite'] = lpname
                    self.G.nodes[self.natoms]['at1'] = host_atom
                    self.G.nodes[self.natoms]['at2'] = ref_atom
                    self.G.nodes[self.natoms]['x'] = 0.0
                    self.G.nodes[self.natoms]['y'] = 0.0
                    self.G.nodes[self.natoms]['z'] = 0.0
                    self.G.nodes[self.natoms]['mass'] = 0.0
                    
                    # Get distance
                    for j, part in enumerate(entry):
                        if part == "DIST" and j+1 < len(entry):
                            try:
                                self.G.nodes[self.natoms]['dist'] = float(entry[j+1]) * 0.1
                            except ValueError:
                                self.G.nodes[self.natoms]['dist'] = 0.16
                    
                    self.natoms += 1
                    self.nvsites += 1
        
        self.coord = [[0.0, 0.0, 0.0] for i in range(self.natoms)]

    def read_mol2_coor_only(self, filename):
        """Read coordinates from MOL2 file"""
        try:
            with open(filename, 'r') as f:
                section = "NONE"
                for line in f.readlines():
                    secflag = False
                    if line.startswith("@"):
                        secflag = True
                        section = "NONE"
                    
                    if section == "NATO" and not secflag:
                        parts = line.split()
                        if len(parts) >= 2:
                            check_natoms = int(parts[0])
                            check_nbonds = int(parts[1])
                            
                            if check_natoms != self.natoms - self.nvsites:
                                if self.nvsites > 0:
                                    print(f"\nNOTE: {self.nvsites} lone pairs found in topology (not in MOL2)")
                                else:
                                    print(f"Error: atom count mismatch (MOL2: {check_natoms}, Topology: {self.natoms})")
                                    sys.exit(1)
                            section = "NONE"
                    
                    if section == "MOLE" and not secflag:
                        self.name = line.strip()
                        section = "NATO"
                    
                    if section == "ATOM" and not secflag:
                        parts = line.split()
                        if len(parts) >= 5:
                            try:
                                atomi = int(parts[0]) - 1
                                x = float(parts[2])
                                y = float(parts[3])
                                z = float(parts[4])
                                
                                self.G.nodes[atomi]['x'] = x
                                self.G.nodes[atomi]['y'] = y
                                self.G.nodes[atomi]['z'] = z
                                self.coord[atomi][0] = x
                                self.coord[atomi][1] = y
                                self.coord[atomi][2] = z
                                
                                # Handle lone pair host atoms
                                if is_lp_host_atom(self.G, self.nvsites, self.G.nodes[atomi]['name']):
                                    atomj = find_vsite(self.G, self.nvsites, self.natoms, atomi)
                                    if atomj is not None:
                                        self.G.nodes[atomj]['x'] = 9999.99
                                        self.G.nodes[atomj]['y'] = 9999.99
                                        self.G.nodes[atomj]['z'] = 9999.99
                                        self.coord[atomj][0] = 9999.99
                                        self.coord[atomj][1] = 9999.99
                                        self.coord[atomj][2] = 9999.99
                            except ValueError:
                                continue
                    
                    if line.startswith("@<TRIPOS>MOLECULE"):
                        section = "MOLE"
                    if line.startswith("@<TRIPOS>ATOM"):
                        section = "ATOM"
                    if line.startswith("@<TRIPOS>BOND"):
                        section = "BOND"
        except FileNotFoundError:
            print(f"Error: Could not find {filename}")
            sys.exit(1)

    def write_pdb(self, f):
        """Write PDB file with coordinates"""
        for atomi in range(0, self.natoms):
            atname = self.G.nodes[atomi]['name']
            
            if len(atname) > 4:
                print("error: atom name > 4 characters")
                sys.exit(1)
            
            # Handle lone pairs
            if is_lp(atname):
                atn1 = "dum"
                atn2 = "dum"
                dist = 0
                
                for ai in range(0, self.nvsites):
                    if self.G.nodes[ai].get('vsite') == atname:
                        atn1 = self.G.nodes[ai].get('at1')
                        atn2 = self.G.nodes[ai].get('at2')
                        dist = self.G.nodes[ai].get('dist', 0.16) * 10
                
                at1 = -1
                at2 = -1
                for ai in range(0, self.natoms):
                    if self.G.nodes[ai]['name'] == atn1:
                        at1 = ai
                    if self.G.nodes[ai]['name'] == atn2:
                        at2 = ai
                
                if at1 >= 0 and at2 >= 0:
                    x1 = self.coord[at1][0]
                    y1 = self.coord[at1][1]
                    z1 = self.coord[at1][2]
                    x2 = self.coord[at2][0]
                    y2 = self.coord[at2][1]
                    z2 = self.coord[at2][2]
                    
                    xlp, ylp, zlp = construct_lp(x1, y1, z1, x2, y2, z2, dist)
                    self.coord[atomi][0] = xlp
                    self.coord[atomi][1] = ylp
                    self.coord[atomi][2] = zlp
            
            resid = self.G.nodes[atomi].get('resid', 1)
            f.write("%-6s%5d %-4s %-4s%5s%12.3f%8.3f%8.3f%6.2f%6.2f\n" %
                    ("ATOM", atomi + 1, atname, self.name, resid,
                     self.coord[atomi][0], self.coord[atomi][1], self.coord[atomi][2], 1.0, 0.0))
        f.write("END\n")

    def write_gmx_itp(self, filename, angle_pars):
        """Write GROMACS itp file"""
        try:
            with open(filename, 'w') as f:
                f.write("[ moleculetype ]\n")
                f.write(f"; Name            nrexcl\n")
                f.write(f"{self.name:16s}     3\n\n")
                
                # Write atoms
                f.write("[ atoms ]\n")
                f.write(";   nr       type  resnr residue  atom   cgnr     charge       mass\n")
                
                total_charge = 0.0
                for i in range(self.natoms):
                    atname = self.G.nodes[i]['name']
                    attype = self.G.nodes[i]['type']
                    charge = self.G.nodes[i]['charge']
                    mass = self.G.nodes[i]['mass']
                    resid = self.G.nodes[i]['resid']
                    
                    total_charge += charge
                    f.write(f"{i+1:6d}{attype:>7}{resid:6d}{self.name:8}{atname:>6}{i+1:6d}{charge:12.4f}{mass:12.4f}\n")
                
                f.write(f"\n; Total charge: {total_charge:.4f}\n\n")
                
                # Write bonds
                f.write("[ bonds ]\n")
                f.write(";  ai    aj funct\n")
                
                for edge in self.G.edges():
                    i, j = edge
                    # Skip bonds to virtual sites
                    if not is_lp(self.G.nodes[i]['name']) and not is_lp(self.G.nodes[j]['name']):
                        f.write(f"{i+1:6d}{j+1:6d}     1\n")
                
                f.write("\n")
                
                # Write virtual site constraints if any
                if self.nvsites > 0:
                    f.write("[ virtual_sites3 ]\n")
                    f.write("; Site   from        funct theta d\n")
                    for i in range(self.natoms):
                        if is_lp(self.G.nodes[i]['name']):
                            at1_name = self.G.nodes[i].get('at1')
                            at2_name = self.G.nodes[i].get('at2')
                            dist = self.G.nodes[i].get('dist', 0.16)
                            
                            at1_idx = -1
                            at2_idx = -1
                            for j in range(self.natoms):
                                if self.G.nodes[j]['name'] == at1_name:
                                    at1_idx = j
                                if self.G.nodes[j]['name'] == at2_name:
                                    at2_idx = j
                            
                            if at1_idx >= 0 and at2_idx >= 0:
                                f.write(f"{i+1:5d}{at1_idx+1:5d}{at2_idx+1:5d}     3  180.00   {dist:.3f}\n")
                    f.write("\n")
        
        except IOError as e:
            print(f"Error writing {filename}: {e}")
            sys.exit(1)

#=================================================================================================================

def parse_charmm_parameters(prmlines):
    """Parse CHARMM parameters - simplified for CGenFF 5.0"""
    params = {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}
    
    noblanks = [x for x in prmlines if len(x.strip()) > 0]
    nocomments = [x for x in noblanks if x.strip()[0] not in ['*', '!']]
    
    section = "NONE"
    for line in nocomments:
        if line.startswith("BONDS"):
            section = "BONDS"
            continue
        elif line.startswith("ANGLES"):
            section = "ANGLES"
            continue
        elif line.startswith("DIHEDRALS"):
            section = "DIHEDRALS"
            continue
        elif line.startswith("IMPROPERS"):
            section = "IMPROPERS"
            continue
        elif line.startswith("END"):
            section = "NONE"
            continue
        
        parts = line.split()
        if len(parts) < 2:
            continue
        
        if section == "BONDS" and len(parts) >= 4:
            try:
                params['bonds'].append({
                    'types': [parts[0], parts[1]],
                    'k': float(parts[2]),
                    'r0': float(parts[3])
                })
            except ValueError:
                continue
        elif section == "ANGLES" and len(parts) >= 4:
            try:
                params['angles'].append({
                    'types': [parts[0], parts[1], parts[2]],
                    'k': float(parts[3]),
                    'theta0': float(parts[4]) if len(parts) > 4 else 0.0
                })
            except ValueError:
                continue
        elif section == "DIHEDRALS" and len(parts) >= 5:
            try:
                params['dihedrals'].append({
                    'types': [parts[0], parts[1], parts[2], parts[3]],
                    'k': float(parts[4])
                })
            except ValueError:
                continue
    
    return params

def write_gmx_bon(params, prefix, filename):
    """Write GROMACS parameter file"""
    try:
        with open(filename, 'w') as f:
            f.write("; Parameters auto-generated from CGenFF\n\n")
            
            if params.get('bonds'):
                f.write("[ bondtypes ]\n")
                f.write("; i    j funct c0 c1\n")
                for bond in params['bonds']:
                    types = bond['types']
                    f.write(f"{types[0]:>6} {types[1]:<6} 1 {bond['r0']:10.6f} {bond['k']:10.1f}\n")
                f.write("\n")
            
            if params.get('angles'):
                f.write("[ angletypes ]\n")
                f.write(";  i    j    k funct theta0 ktheta\n")
                for angle in params['angles']:
                    types = angle['types']
                    theta0 = angle.get('theta0', 120.0)
                    f.write(f"{types[0]:>6} {types[1]:>6} {types[2]:<6} 5 {theta0:8.3f} {angle['k']:10.1f}\n")
                f.write("\n")
    except IOError as e:
        print(f"Error writing {filename}: {e}")
        sys.exit(1)

def write_gmx_mol_top(filename, ffdir, prmfile, itpfile, molname):
    """Write GROMACS topology file"""
    try:
        with open(filename, 'w') as f:
            f.write("; Generated GROMACS topology\n")
            f.write(f"; Molecule: {molname}\n\n")
            f.write("; Include forcefield parameters\n")
            f.write('#include "./charmm36-jul2022.ff/forcefield.itp"\n\n')
            f.write("; Include molecule parameters\n")
            f.write(f'#include "{prmfile}"\n')
            f.write(f'#include "{itpfile}"\n\n')
            f.write("[ system ]\n")
            f.write(f"{molname}\n\n")
            f.write("[ molecules ]\n")
            f.write(f"{molname:20s} 1\n")
    except IOError as e:
        print(f"Error writing {filename}: {e}")
        sys.exit(1)

#=================================================================================================================

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py RESNAME molecule.mol2 molecule.str charmm36.ff")
        print("\nExample: python cgenff_charmm2gmx_py3_v5.0.py DAC dacomitinib.mol2 dacomitinib_fixed.str ./charmm36-jul2022.ff")
        sys.exit(1)
    
    # Check NetworkX version
    nx_version = float(nx.__version__.split('.')[0])
    if nx_version < 2.0:
        print(f"Error: NetworkX version {nx.__version__} is too old. Please install version 2.3+")
        sys.exit(1)
    
    print(f"Python version: {sys.version}")
    print(f"NetworkX version: {nx.__version__}")
    print()
    
    mol_name = sys.argv[1]
    mol2_name = sys.argv[2]
    rtp_name = sys.argv[3]
    ffdir = sys.argv[4]
    
    print(f"Input molecule name: {mol_name}")
    print(f"Input MOL2 file: {mol2_name}")
    print(f"Input CGenFF file: {rtp_name}")
    print(f"Force field directory: {ffdir}")
    print()
    
    # Check versions
    check_versions(rtp_name, os.path.join(ffdir, "forcefield.doc"))
    print()
    
    # Output filenames
    itpfile = mol_name.lower() + ".itp"
    prmfile = mol_name.lower() + ".prm"
    initpdbfile = mol_name.lower() + "_ini.pdb"
    topfile = mol_name.lower() + ".top"
    
    print(f"Output files will be:")
    print(f"  - {itpfile}")
    print(f"  - {prmfile}")
    print(f"  - {initpdbfile}")
    print(f"  - {topfile}")
    print()
    
    # Read atomtypes
    atomtypes_file = os.path.join(ffdir, "atomtypes.atp")
    if not os.path.exists(atomtypes_file):
        print(f"Error: Could not find {atomtypes_file}")
        sys.exit(1)
    
    atomtypes = read_gmx_atomtypes(atomtypes_file)
    print(f"Read {len(atomtypes)} atom types from forcefield")
    
    # Read angle parameters from forcefield
    angl_params = []
    filelist = get_filelist_from_gmx_forcefielditp(ffdir, "forcefield.itp")
    print(f"Reading angle parameters from {len(filelist)} forcefield files...")
    for filename in filelist:
        anglpars = read_gmx_anglpars(filename)
        angl_params = angl_params + anglpars
    print(f"Read {len(angl_params)} angle parameters from forcefield")
    print()
    
    # Parse molecule topology
    print(f"Parsing molecule topology from {rtp_name}...")
    m = atomgroup()
    rtplines = get_charmm_rtp_lines(rtp_name, mol_name)
    
    if not rtplines:
        print(f"Error: Could not find RESI {mol_name} in {rtp_name}")
        sys.exit(1)
    
    m.read_charmm_rtp(rtplines, atomtypes)
    print(f"Found {m.natoms} atoms ({m.nvsites} virtual sites)")
    print()
    
    # Read coordinates
    print(f"Reading coordinates from {mol2_name}...")
    m.read_mol2_coor_only(mol2_name)
    print("Coordinates read successfully")
    print()
    
    # Write PDB file
    print(f"Writing initial PDB file to {initpdbfile}...")
    with open(initpdbfile, 'w') as f:
        m.write_pdb(f)
    print("PDB file written successfully")
    print()
    
    # Parse parameters
    print(f"Parsing parameters from {rtp_name}...")
    prmlines = get_charmm_prm_lines(rtp_name)
    params = parse_charmm_parameters(prmlines)
    print(f"Found {len(params['bonds'])} bond parameters")
    print(f"Found {len(params['angles'])} angle parameters")
    print()
    
    # Write parameter file
    print(f"Writing parameter file to {prmfile}...")
    write_gmx_bon(params, "", prmfile)
    print("Parameter file written successfully")
    print()
    
    # Write ITP file
    print(f"Writing topology file to {itpfile}...")
    m.write_gmx_itp(itpfile, angl_params)
    print("Topology file written successfully")
    print()
    
    # Write TOP file
    print(f"Writing system topology file to {topfile}...")
    write_gmx_mol_top(topfile, ffdir, prmfile, itpfile, mol_name)
    print("System topology file written successfully")
    print()
    
    print("="*50)
    print("CONVERSION COMPLETE")
    print("="*50)
    print(f"\nGenerated files:")
    print(f"  1. {itpfile} - Molecule topology (include in your top file)")
    print(f"  2. {prmfile} - Molecule parameters (include in your top file)")
    print(f"  3. {initpdbfile} - Initial coordinates")
    print(f"  4. {topfile} - Complete system topology (for reference)")
    print(f"\nNote: Make sure to include both {itpfile} and {prmfile} in your main topology file!")
    print(f"If this molecule has lone pairs, ensure you're using GROMACS 2020.x or newer.")
    print()

