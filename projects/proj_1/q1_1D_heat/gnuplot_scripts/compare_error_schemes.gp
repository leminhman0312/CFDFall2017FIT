# gnuplot_scripts/compare_error_schemes.gp

set terminal pngcairo size 1000,700 enhanced font "Arial,18"
set output outpng

set datafile separator "\t"
set datafile commentschars "#"

set xlabel "X [ft]"
set ylabel "Error = T_scheme - T_exact [degree F]"
set title sprintf("Error vs Exact, dt = %.2f hr, t = %s hr", dt, tlabel)

set grid
set key outside right center
set xrange [0:1]

plot \
  0 with lines lw 6 dt 2 lc rgb "black" title "Exact error (0)", \
  sprintf("data/error_ftcs_explicit_%s.txt", tag) index idx using 1:2 \
    with linespoints lw 1.2 pt 6  ps 1.3 lc rgb "red"    title "FTCS explicit", \
  sprintf("data/error_ftcs_implicit_%s.txt", tag) index idx using 1:2 \
    with linespoints lw 1.2 pt 4  ps 1.3 lc rgb "blue"   title "FTCS implicit", \
    sprintf("data/error_dufort_%s.txt", tag) index idx using 1:2 \
    with linespoints lw 1.2 pt 8 ps 1.3 lc rgb "#006400" title "Dufort Frankel", \
  sprintf("data/error_cn_%s.txt", tag) index idx using 1:2 \
    with linespoints lw 1.2 pt 10 ps 1.3 lc rgb "purple" title "Crank Nicolson"
