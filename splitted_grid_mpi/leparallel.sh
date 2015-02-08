#!/bin/bash

grid=200
time=100
n=24
step=3

for ((process=1;process<=$n;process+=$step))
do
	echo $process
	echo -n $process >> le2d_mpi.txt
	echo -n " " >> le2d_mpi.txt
	mpirun -np $process ./le2d_mpi $grid $grid $time >> le2d_mpi.txt
done

echo "OK"
