#!/bin/sh
#sed_anchor01
#SBATCH --output=B13N2_mp-534-very_fine.sout
#SBATCH --job-name=B13N2_mp-534-very_fine
#SBATCH --nodes=1
#SBATCH --partition=C16M32
###SBATCH --exclude=node[01-03,05-07,20-23,29,30]
#SBATCH --ntasks-per-node=16

export LD_LIBRARY_PATH=/opt/mpich-3.4.2/lib:/opt/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
export PATH=/opt/mpich-3.4.2/bin:$PATH
#sed_anchor02
mpiexec /opt/QEGCC_MPICH3.4.2/bin/pw.x -in B13N2_mp-534.in



#squeue       # qstat
#scancel  1   # Delete JOB NO. 1 
#scontrol show job 1  # Check job NO.1 detail



################################################## DEBUG
#yum install dos2unix  #######   install this when you got as follows error messages
#sbatch: error: Batch script contains DOS line breaks (\r\n)
#sbatch: error: instead of expected UNIX line breaks (\n).
#dos2unix filename 
