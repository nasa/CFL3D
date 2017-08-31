#!/bin/bash
#SBATCH --job-name=test_cfl3d         # Job name
#SBATCH --partition=bigdata           # running on bigdata partition
#SBATCH --ntasks=24                   # Run on all CPUs of a node
#SBATCH --nodes=1                     # Use one node
#SBATCH --mem=2500mb                  # Memory limit, 
#SBATCH --time=00:25:00               # Time limit hrs:min:sec
#SBATCH --output=test_cfl3d.out       # Standard output and error log

echo "-----------------Environment------------------"
pwd
hostname
env
echo "---------------------Start--------------------"
date
echo "-----------------CGNS test--------------------"
cd Delta_cgns
rm cfl3d.error -f
rm cfl3d.out -f
plot3dg_to_cgns < plot3dg_to_cgns.inp
cfl3d_seq < delta_cgns.inp
cat cfl3d.error
cd ..
echo "------------------MPI test--------------------"
cd Axibump
rm cfl3d.error -f
rm cfl3d.out -f
splitter < split.inp_3blk
srun -n 4 cfl3d_mpi < bumpperiodic.inp_3blk
cat cfl3d.error
echo "---------------------End----------------------"
date
