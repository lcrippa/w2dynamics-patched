#!/bin/bash
#PBS -N EDcheck_CIRO
#PBS -V
#PBS -l nodes=2:ppn=12
#PBS -l walltime=168:00:00

echo "Starting script..."
BIN=$PBS_O_HOME/opt/bin

cd $PBS_O_WORKDIR

#cat $PBS_NODEFILE > hostlist.$PBS_JOBID
NPROCS=`wc -l < $PBS_NODEFILE`

trap '' USR1 USR2
mpirun -v -machinefile $PBS_NODEFILE -np $NPROCS $BIN/python $BIN/DMFT.py >ctqmc-$PBS_JOBID.log 2>&1

echo "Script complete."
