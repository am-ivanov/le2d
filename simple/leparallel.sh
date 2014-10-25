for ((process=1;process <= 24;process+=2))
do
	for ((i=0;i<5;i++))
	do
		./simple_gcc 4000 4000 100 $process >> simple_gcc.txt
		./simple 4000 4000 100 $process >> simple.txt
	done
	echo "" >> simple_gcc.txt
	echo "" >> simple.txt
done
echo "OK"
