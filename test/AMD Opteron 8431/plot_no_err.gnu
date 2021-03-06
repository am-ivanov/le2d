set term png size 1024,768
set output "speedup_no_err.png"
set style data lines
set grid xtics ytics mxtics mytics
set title "amd 48 cores"
set xlabel "Threads"
set ylabel "Speedup"

plot "spl_gr_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "splitted grid", \
"simple_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "simple", \
"cf_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "cache friendly", \
"pthread_gcc_toplot.txt" using ($0*4 + 1):($1/$3):(sqrt((($2)**2) + (($4)**2))) ti "pthreads"
