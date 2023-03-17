#!/bin/bash

# Set the name of the executable
OBJ=a.out

# A command to execute a parallel program with MPI
MPIRUN=mpirun

# number of processes to be spawned
NUM=`grep "N_beta" input.ini | awk '{print$2}'`
echo ""
echo "MSG from mjob.sh: We will run" $NUM "(N_beta) running processes."
echo ""

# Run the executable file
#
# to default to the number of hardware threads instead of the number of processor cores
$MPIRUN -n $NUM --use-hwthread-cpus  ./$OBJ
#
#  to ignore the number of available slots when deciding the number of processes to launch
#$MPIRUN -n $NUM --oversubscribe ./$OBJ
#
#$MPIRUN -n $NUM ./$OBJ

# Delete the executable file
#rm $OBJ
