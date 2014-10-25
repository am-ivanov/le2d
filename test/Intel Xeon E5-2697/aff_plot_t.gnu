set term png
set terminal png size 1024,768
set output "aff_time.png"
set style data lines
set grid xtics ytics mxtics mytics
set title "intel 24 cores + affinity"
set xlabel "Threads"
set ylabel "Time"
set xrange [0:25]

plot "aff_spl_gr_toplot.txt" using ($0*2 + 1):($3):($4) ti "splitted grid icc no vect", \
"aff_spl_gr_gcc_toplot.txt" using ($0*2 + 1):($3):($4) ti "splitted grid gcc no vect", \
"aff_spl_gr_av_toplot.txt" using ($0*2 + 1):($3):($4) ti "splitted grid icc vect", \
"aff_cf_toplot.txt" using ($0*2 + 1):($3):($4) ti "cache friendly icc no vect", \
"aff_cf_gcc_toplot.txt" using ($0*2 + 1):($3):($4) ti "cache friendly gcc no vect", \
"aff_cf_av_toplot.txt" using ($0*2 + 1):($3):($4) ti "cache friendly icc vect", \
"aff_pthread_toplot.txt" using ($0*2 + 1):($3):($4) ti "pthread icc no vect", \
"aff_pthread_gcc_toplot.txt" using ($0*2 + 1):($3):($4) ti "pthread gcc no vect", \
"aff_pthread_av_toplot.txt" using ($0*2 + 1):($3):($4) ti "pthread icc vect", \
"aff_simple_gcc_toplot.txt" using ($0*2 + 1):($3):($4) ti "simple gcc no vect", \
"aff_simple_toplot.txt" using ($0*2 + 1):($3):($4) ti "simple icc no vect"
