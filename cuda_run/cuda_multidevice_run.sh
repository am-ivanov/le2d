#!/bin/bash

x=4096
y=4096
t=6500

for file in cuda_le2d_f cuda_le2d_f_fm cuda_le2d_d 
do
	for ((vers=5;vers<=6;++vers))
	do
		CUDA_VISIBLE_DEVICES=0 cuda$vers/$file $x $y $t >> result_multidevice.txt
		#CUDA_VISIBLE_DEVICES=0,1 cuda$vers/$file $x $y $t >> result_multidevice.txt
		#CUDA_VISIBLE_DEVICES=0,1,2 cuda$vers/$file $x $y $t >> result_multidevice.txt
		#CUDA_VISIBLE_DEVICES=0,1,2,3 cuda$vers/$file $x $y $t >> result_multidevice.txt
		#CUDA_VISIBLE_DEVICES=0,1,2,3,4 cuda$vers/$file $x $y $t >> result_multidevice.txt
	done
done

echo OK
