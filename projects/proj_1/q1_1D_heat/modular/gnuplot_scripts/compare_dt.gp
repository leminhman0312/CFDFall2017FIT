# usage:
# gnuplot -e "scheme='ftcs_implicit'; tag1='001'; dt1='0.01'; tag2='005'; dt2='0.05'; idx=1; tlabel='0.1'" compare_dt.gp

set terminal pngcairo size 1000,700 enhanced font "Arial,18"
set output sprintf("compare/compare_dt_%s_t%s.png", scheme, tlabel)

set datafile separator "\t"
set datafile commentschars "#"

set xlabel "X [ft]"
set ylabel "T [degree F]"

sname = scheme
if (scheme eq "ftcs_explicit") sname = "FTCS explicit"
if (scheme eq "ftcs_implicit") sname = "FTCS implicit"
if (scheme eq "dufort")        sname = "Dufort-Frankel"
if (scheme eq "cn")            sname = "Crank-Nicolson"

set title sprintf("%s: compare dt at t = %s hr", sname, tlabel)

set grid
set key outside right center
set xrange [0:1]
set yrange [100:300]

f1 = sprintf("data/%s_%s.txt", scheme, tag1)
f2 = sprintf("data/%s_%s.txt", scheme, tag2)

plot \
  f1 index idx using 1:2 with linespoints lw 2 pt 7  ps 1.0 title sprintf("dt = %s hr", dt1), \
  f2 index idx using 1:2 with linespoints lw 2 pt 5  ps 1.0 title sprintf("dt = %s hr", dt2)