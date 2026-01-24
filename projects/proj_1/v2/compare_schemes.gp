# usage:
# gnuplot -e "tag='001'; dt='0.01'; idx=1; tlabel='0.1'" compare_schemes.gp

set terminal pngcairo size 1000,700 enhanced font "Arial,18"
set output sprintf("compare_schemes_dt%s_t%s.png", tag, tlabel)

set datafile separator "\t"
set datafile commentschars "#"

set xlabel "X [ft]"
set ylabel "T [degree F]"
set title sprintf("Compare schemes at dt = %s hr, t = %s hr", dt, tlabel)

set grid
set key outside right center
set xrange [0:1]
set yrange [100:300]

plot \
  sprintf("ftcs_explicit_%s.txt", tag) index idx using 1:2 with linespoints lw 2 pt 7  ps 1.0 title "FTCS explicit", \
  sprintf("dufort_%s.txt",        tag) index idx using 1:2 with linespoints lw 2 pt 5  ps 1.0 title "Dufort-Frankel", \
  sprintf("ftcs_implicit_%s.txt", tag) index idx using 1:2 with linespoints lw 2 pt 9  ps 1.0 title "FTCS implicit", \
  sprintf("cn_%s.txt",            tag) index idx using 1:2 with linespoints lw 2 pt 11 ps 1.0 title "Crank-Nicolson"
