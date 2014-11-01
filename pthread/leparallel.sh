n=24
step=2
grid=4000
time=100

unset GOMP_CPU_AFFINITY

for ((process=1;process<=$n;process+=$step))
do
	for ((i=0;i<5;i++))
	do
		./pthread_gcc_f $grid $grid $time $process >> pthread_gcc_f.txt
		./pthread_icc_f $grid $grid $time $process >> pthread_icc_f.txt
		./pthread_gcc_d $grid $grid $time $process >> pthread_gcc_d.txt
		./pthread_icc_d $grid $grid $time $process >> pthread_icc_d.txt
		./pthread_icc_sse_f $grid $grid $time $process >> pthread_icc_sse_f.txt
		./pthread_icc_sse_d $grid $grid $time $process >> pthread_icc_sse_d.txt
		./pthread_icc_avx_f $grid $grid $time $process >> pthread_icc_avx_f.txt
		./pthread_icc_avx_d $grid $grid $time $process >> pthread_icc_avx_d.txt
	done
	echo "" >> pthread_gcc_f.txt
	echo "" >> pthread_icc_f.txt
	echo "" >> pthread_gcc_d.txt
	echo "" >> pthread_icc_d.txt
	echo "" >> pthread_icc_sse_f.txt
	echo "" >> pthread_icc_sse_d.txt
	echo "" >> pthread_icc_avx_f.txt
	echo "" >> pthread_icc_avx_d.txt
done

echo "OK"