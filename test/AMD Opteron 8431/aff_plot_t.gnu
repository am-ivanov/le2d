set term png
set terminal png size 1024,768 
set output "aff_time.png"
set style data lines
set grid xtics ytics mxtics mytics
set title "amd 48 cores + affinity"
set xlabel "Threads"
set ylabel "Time"

plot "aff_spl_gr_gcc_toplot.txt" using ($0*4 + 1):($3):($4) ti "splitted grid", \
"aff_simple_gcc_toplot.txt" using ($0*4 + 1):($3):($4) ti "simple", \
"aff_cf_gcc_toplot.txt" using ($0*4 + 1):($3):($4) ti "cahce friendly", \
"aff_pthread_gcc_toplot.txt" using ($0*4 + 1):($3):($4) ti "pthreads"
