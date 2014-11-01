target=gcc_compiler

echo comiling binaries
make -f ../simple/Makefile -C ../simple $target
make -f ../cache_friendly/Makefile -C ../cache_friendly $target
make -f ../pthread/Makefile -C ../pthread $target
make -f ../splitted_grid/Makefile -C ../splitted_grid $target

echo copying in current directory
cp ../simple/simple_* ./
cp ../cache_friendly/cf_* ./
cp ../pthread/pthread_* ./
cp ../splitted_grid/spl_gr_* ./

echo cleaning source directories
make -f ../simple/Makefile -C ../simple clean
make -f ../cache_friendly/Makefile -C ../cache_friendly clean
make -f ../pthread/Makefile -C ../pthread clean
make -f ../splitted_grid/Makefile -C ../splitted_grid clean

echo copying scripts
cp ../simple/leparallel.sh ./leparallel_simple.sh
cp ../cache_friendly/leparallel.sh ./leparallel_cf.sh
cp ../pthread/leparallel.sh ./leparallel_pthread.sh
cp ../splitted_grid/leparallel.sh ./leparallel_spl_gr.sh

echo test start
echo simple version
./leparallel_simple.sh
echo cache friendly version
./leparallel_cf.sh
echo pthread version
./leparallel_pthread.sh
echo splitted grid version
./leparallel_spl_gr.sh

if (($target == gcc_compiler)) 
	then
		rm *_icc_*
fi

echo "OK"