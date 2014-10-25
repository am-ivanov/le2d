set term png size 1024,768
set output "aff_speedup_no_err.png"
set style data lines
set grid xtics ytics mxtics mytics
set title "amd 48 cores + affinity"
set xlabel "Threads"
set ylabel "Speedup"

plot "aff_spl_gr_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "splitted grid", \
"aff_simple_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "simple", \
"aff_cf_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "cache friendly", \
"aff_pthread_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "pthreads"
