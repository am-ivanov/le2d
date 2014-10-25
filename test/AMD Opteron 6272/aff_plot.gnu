set term png
set terminal png size 1024,768
set output "aff_speedup.png"
set style data yerrorbars
set grid
set title "amd 64 cores + affinity"
set xlabel "Threads"
set ylabel "Speedup"

plot "aff_spl_gr_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "splitted grid",\
"aff_simple_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "simple",\
"aff_cf_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "cache friendly",\
"aff_pthread_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "pthreads"
