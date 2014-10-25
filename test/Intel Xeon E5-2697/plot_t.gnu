set term png
set terminal png size 1024,768
set output "time.png"
set style data lines
set grid xtics ytics mxtics mytics
set title "intel 24 cores"
set xlabel "Threads"
set ylabel "Time"
set xrange [0:25]

plot "spl_gr_toplot.txt" using ($0*2 + 1):($3):($4) ti "splitted grid icc no vect", \
"spl_gr_gcc_toplot.txt" using ($0*2 + 1):($3):($4) ti "splitted grid gcc no vect", \
"spl_gr_av_toplot.txt" using ($0*2 + 1):($3):($4) ti "splitted grid icc vect", \
"cf_toplot.txt" using ($0*2 + 1):($3):($4) ti "cache friendly icc no vect", \
"cf_gcc_toplot.txt" using ($0*2 + 1):($3):($4) ti "cache friendly gcc no vect", \
"cf_av_toplot.txt" using ($0*2 + 1):($3):($4) ti "cache friendly icc vect", \
"pthread_toplot.txt" using ($0*2 + 1):($3):($4) ti "pthread icc no vect", \
"pthread_gcc_toplot.txt" using ($0*2 + 1):($3):($4) ti "pthread gcc no vect", \
"pthread_av_toplot.txt" using ($0*2 + 1):($3):($4) ti "pthread icc vect", \
"simple_gcc_toplot.txt" using ($0*2 + 1):($3):($4) ti "simple gcc no vect", \
"simple_toplot.txt" using ($0*2 + 1):($3):($4) ti "simple icc no vect"