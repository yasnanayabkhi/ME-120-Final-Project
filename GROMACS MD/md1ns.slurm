#!/bin/bash

#================================================================================
# md1ns.slurm
#
# GROMACS 1 ns Production MD for FtsZ–Dacomtinib Complex
#================================================================================

#SBATCH --job-name=md
#SBATCH --output=logs/md.out
#SBATCH --error=logs/md.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yasna_nayabkhil@berkeley.edu
#SBATCH --mem=3500M
#SBATCH --time=2-00:00:00
#SBATCH --exclusive

#=============================================================================
# Load modules

module purge
module load ohpc
module load spack/1.0.2
module load gnu15
module load openmpi5/5.0.7
module load gromacs/2024.5-gcc-11.5.0

#=============================================================================
# MPI settings (disable infiniband since cluster uses ethernet)

export OMPI_MCA_btl="^openib"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#=============================================================================
# Job info

echo "======================================"
echo "Starting 1 ns Production MD"
echo "Node: $(hostname)"
echo "MPI ranks: $SLURM_NTASKS"
echo "Start time: $(date)"
echo "======================================"

#=============================================================================
# Generate MD run file

gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1

#=============================================================================
# Run MD simulation

mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm md -v

#=============================================================================

echo "======================================"
echo "1 ns MD production finished"
echo "End time: $(date)"
echo "======================================"