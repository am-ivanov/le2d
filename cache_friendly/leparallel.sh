for ((process=1;process <= 24;process+=2))
do
	for ((i=0;i<5;i++))
	do
		./cf_gcc 4000 4000 100 $process >> cf_gcc.txt
		./cf 4000 4000 100 $process >> cf.txt
		./cf_av 4000 4000 100 $process >> cf_av.txt
	done
	echo "" >> cf_gcc.txt
	echo "" >> cf.txt
	echo "" >> cf_av.txt
done
echo "OK"
