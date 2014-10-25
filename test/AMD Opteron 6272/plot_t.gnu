set term png
set terminal png size 1024,768 
set output "time.png"
set style data lines
set grid xtics ytics mxtics mytics
set title "amd 64 cores"
set xlabel "Threads"
set ylabel "Time"

plot "spl_gr_gcc_toplot.txt" using ($0*4 + 1):($3):($4) ti "splitted grid", \
"simple_gcc_toplot.txt" using ($0*4 + 1):($3):($4) ti "simple", \
"cf_gcc_toplot.txt" using ($0*4 + 1):($3):($4) ti "cahce friendly", \
"pthread_gcc_toplot.txt" using ($0*4 + 1):($3):($4) ti "pthreads"
