for ((process=1;process <= 24;process+=2))
do
	for ((i=0;i<5;i++))
	do
		./spl_gr_gcc 4000 4000 100 $process >> spl_gr_gcc.txt
		./spl_gr 4000 4000 100 $process >> spl_gr.txt
		./spl_gr_av 4000 4000 100 $process >> spl_gr_av.txt
	done
	echo "" >> spl_gr_gcc.txt
	echo "" >> spl_gr.txt
	echo "" >> spl_gr_av.txt
done
echo "OK"
