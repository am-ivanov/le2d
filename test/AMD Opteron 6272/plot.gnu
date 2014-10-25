set term png
set terminal png size 1024,768
set output "speedup.png"
set style data yerrorbars
set grid
set title "amd 64 cores"
set xlabel "Threads"
set ylabel "Speedup"

plot "spl_gr_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "splitted frid",\
"simple_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "simple",\
"cf_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "cache friendly",\
"pthread_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "pthreads"
