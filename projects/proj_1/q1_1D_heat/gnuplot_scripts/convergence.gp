set terminal pngcairo size 1000,700 enhanced font "Arial,18"
set output outpng

set datafile separator "\t"
set datafile commentschars "#"

set grid
set key outside right center

set xlabel "dt [hr]"
set ylabel "Error norm"
set title sprintf("Convergence vs dt at t = %s hr", tlabel)

set logscale x
set logscale y

plot \
  infile using 1:2 with linespoints lw 2 pt 6 lc rgb "#d62728" title "FTCS explicit L2", \
  infile using 1:3 with linespoints lw 2 pt 4 lc rgb "#1f77b4" title "FTCS implicit L2", \
  infile using 1:4 with linespoints lw 2 pt 8 lc rgb "#2ca02c" title "Dufort Frankel L2", \
  infile using 1:5 with linespoints lw 2 pt 10 lc rgb "#9467bd" title "Crank Nicolson L2"
