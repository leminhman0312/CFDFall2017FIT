set terminal pngcairo size 1800,900 enhanced font "Arial,22"
set output outpng

set xlabel "X [ft]"
set ylabel "Y [ft]"

set title sprintf("Solving {/Symbol \\266}U/{/Symbol \\266}t = {/Symbol a} ( {/Symbol \\266}^2U/{/Symbol \\266}x^2 + {/Symbol \\266}^2U/{/Symbol \\266}y^2 )  at t = %s hr", tlabel)

set size ratio -1
set grid

set xrange [0:xmax]
set yrange [0:ymax]

set pm3d map
set view map
set palette rgbformulae 33,13,10
set colorbox
set cblabel "Temperature"

set contour base
set cntrparam levels 30
unset surface
set key off

splot infile using 1:2:3 with pm3d, \
      infile using 1:2:3 with lines lc rgb "black" lw 0.3
