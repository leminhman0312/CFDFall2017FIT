# gnuplot_scripts/plot_contour_2d_matlab_like.gp

set terminal pngcairo size 1800,800 enhanced font "Arial,22"
set output outpng

set xlabel "X [ft]"

# Make Y label horizontal
set ylabel "Y [ft]" rotate by 0 offset -3,0

# Enhanced text title, pngcairo safe
set title "Solving {/Symbol \266}U/{/Symbol \266}t = {/Symbol a} ( {/Symbol \266}^2U/{/Symbol \266}x^2 + {/Symbol \266}^2U/{/Symbol \266}y^2 ) at t = 0.5 hr"

set size ratio -1
set tics out
set grid back

set view map
set pm3d map

set cntrparam levels 20
set isosamples 200,200

set palette rgbformulae 33,13,10

set cbrange [0:200]
set cbtics (0, 50, 100, 150, 200)

# Make colorbar label horizontal
set cblabel "Temperature" rotate by 0 offset 5,0

unset key

splot datafile using 1:2:3 notitle
