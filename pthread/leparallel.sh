for ((process=1;process <= 24;process+=2))
do
	for ((i=0;i<5;i++))
	do
		./pthread_gcc 4000 4000 100 $process >> pthread_gcc.txt
		./pthread 4000 4000 100 $process >> pthread.txt
		./pthread_av 4000 4000 100 $process >> pthread_av.txt
	done
	echo "" >> pthread_gcc.txt
	echo "" >> pthread.txt
	echo "" >> pthread_av.txt
done
echo "OK"
