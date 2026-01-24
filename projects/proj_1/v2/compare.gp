set terminal pngcairo size 1000,700 enhanced font "Arial,18"
set output "compare.png"

set datafile separator "\t"
set datafile commentschars "#"

set xlabel "X [ft]"
set ylabel "T [degree F]"
set title "{/Symbol \\266}T/{/Symbol \\266}t = {/Symbol a} [{/Symbol \\266}^2T/{/Symbol \\266}x^2] , Compare schemes, dt = 0.01 hr"

set grid
set key outside right center
set xrange [0:1]
set yrange [100:300]



plot \
  "ftcs_explicit.txt"   index 1 using 1:2 with linespoints lw 2 pt 7 ps 1.0 title "FTCS explicit", \
  "dufort.txt" index 1 using 1:2 with linespoints lw 2 pt 5 ps 1.0 title "DuFort Frankel", \
  "cn.txt" index 1 using 1:2 with linespoints lw 2 pt 5 ps 1.0 title "Crank-Nicholson", \
  "ftcs_implicit.txt" index 1 using 1:2 with linespoints lw 2 pt 5 ps 1.0 title "FTCS implicit"