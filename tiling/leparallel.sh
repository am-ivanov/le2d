#!/bin/bash

n=24
step=2
grid=4000
time=100

unset GOMP_CPU_AFFINITY

for ((process=1;process<=$n;process+=$step))
do
	for ((i=0;i<5;i++))
	do
		./simple_gcc_f $grid $grid $time $process >> simple_gcc_f.txt
		./simple_icc_f $grid $grid $time $process >> simple_icc_f.txt
		./simple_gcc_d $grid $grid $time $process >> simple_gcc_d.txt
		./simple_icc_d $grid $grid $time $process >> simple_icc_d.txt
	done
	echo "" >> simple_gcc_f.txt
	echo "" >> simple_icc_f.txt
	echo "" >> simple_gcc_d.txt
	echo "" >> simple_icc_d.txt
done

export GOMP_CPU_AFFINITY="0-$n"

for ((process=1;process<=$n;process+=$step))
do
	for ((i=0;i<5;i++))
	do
		./simple_gcc_f $grid $grid $time $process >> simple_gcc_f_aff.txt
		./simple_icc_f $grid $grid $time $process >> simple_icc_f_aff.txt
		./simple_gcc_d $grid $grid $time $process >> simple_gcc_d_aff.txt
		./simple_icc_d $grid $grid $time $process >> simple_icc_d_aff.txt
	done
	echo "" >> simple_gcc_f_aff.txt
	echo "" >> simple_icc_f_aff.txt
	echo "" >> simple_gcc_d_aff.txt
	echo "" >> simple_icc_d_aff.txt
done

echo "OK"
