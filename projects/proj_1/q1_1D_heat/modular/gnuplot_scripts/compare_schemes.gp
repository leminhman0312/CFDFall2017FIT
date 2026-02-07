set terminal pngcairo size 1000,700 enhanced font "Arial,18"
set datafile separator "\t"
set datafile commentschars "#"

set xlabel "X [ft]"
set ylabel "T [degree F]"
set grid
set key outside right center
set xrange [0:1]
set yrange [100:300]

set output sprintf("compare/compare_schemes_dt%s_t%s.png", tag, tlabel)
set title sprintf("Schemes vs Exact, dt = %s hr, t = %s hr", dt, tlabel)

# exact: bigger black dots only
# schemes: thinner lines, larger point markers, spaced point types
plot \
  sprintf("data/exact_%s.txt", tag) index idx using 1:2 \
    with points pt 7 ps 3.0 lc rgb "black" title sprintf("Exact, t=%s", tlabel), \
  sprintf("data/ftcs_explicit_%s.txt", tag) index idx using 1:2 \
    with linespoints lw 1.2 pt 6 ps 1.6 lc rgb "#1f77b4" title "FTCS explicit", \
  sprintf("data/ftcs_implicit_%s.txt", tag) index idx using 1:2 \
    with linespoints lw 1.2 pt 4 ps 1.6 lc rgb "#d62728" title "FTCS implicit", \
  sprintf("data/dufort_%s.txt", tag) index idx using 1:2 \
    with linespoints lw 1.2 pt 8 ps 1.6 lc rgb "#2ca02c" title "Dufort Frankel", \
  sprintf("data/cn_%s.txt", tag) index idx using 1:2 \
    with linespoints lw 1.2 pt 10 ps 1.6 lc rgb "#9467bd" title "Crank Nicolson"
