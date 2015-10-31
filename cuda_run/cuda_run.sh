#!/bin/bash

x=4096
y=4096
t=6500
gpus=1

for ((vers=1;vers<=4;++vers))
do
	for ((gpu=0;gpu!=gpus;++gpu))
	do
		for file in cuda_le2d_f cuda_le2d_f_fm cuda_le2d_d 
		do
			cuda$vers/$file $x $y $t $gpu >> result.txt
		done
	done
done

echo OK
