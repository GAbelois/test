#!/bin/bash
### Project code
#PBS -A UCLA0034
### Job Name
#PBS -N test
### Merge output and error files
#PBS -j oe
#PBS -q economy

### Select X nodes with Y CPUs and Y mpiprocs each for a total of X*Y MPI processes
#PBS -l select=5:ncpus=32:mpiprocs=32
#PBS -l walltime=12:00:00
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M yanchao@ucla.edu
#PBS -o log.out

export MPI_SHEPHERD=true
export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR
#export ../output/lesgo.log

### Run the executable
module load impi
mpirun -np 160 ./test