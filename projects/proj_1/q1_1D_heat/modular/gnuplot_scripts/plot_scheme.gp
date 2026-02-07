# usage:
# gnuplot -e "scheme='ftcs_explicit'; tag='001'; dt='0.01'; outpng='ftcs_explicit_001.png'" plot_scheme.gp

set terminal pngcairo size 1000,700 enhanced font "Arial,18"
set output outpng

set datafile separator "\t"
set datafile commentschars "#"

set xlabel "X [ft]"
set ylabel "T [degree F]"
# scheme is used for the filename, sname is used for the title
sname = scheme
if (scheme eq "ftcs_explicit") sname = "FTCS explicit"
if (scheme eq "ftcs_implicit") sname = "FTCS implicit"
if (scheme eq "dufort")        sname = "Dufort-Frankel"
if (scheme eq "cn")            sname = "Crank-Nicolson"

set title "{/Symbol \\266}T/{/Symbol \\266}t = {/Symbol a} [{/Symbol \\266}^2T/{/Symbol \\266}x^2] ,  ".sname.", dt = ".dt." hr"

set grid
set key outside right center
set xrange [0:1]
set yrange [100:300]

fname = sprintf("data/%s_%s.txt", scheme, tag)

plot \
  fname index 1 using 1:2 with linespoints lw 2 pt 7  ps 1.0 title "t = 0.1 hr", \
  fname index 2 using 1:2 with linespoints lw 2 pt 5  ps 1.0 title "t = 0.2 hr", \
  fname index 3 using 1:2 with linespoints lw 2 pt 9  ps 1.0 title "t = 0.3 hr", \
  fname index 4 using 1:2 with linespoints lw 2 pt 11 ps 1.0 title "t = 0.4 hr"