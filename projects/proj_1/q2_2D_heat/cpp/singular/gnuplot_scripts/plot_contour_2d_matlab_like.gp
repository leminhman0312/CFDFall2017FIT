# --------------------------------------------------
# 2D Heat contour plot (MATLAB-like)
# --------------------------------------------------

set terminal pngcairo size 1800,800 enhanced font "Arial,22"
set output outpng

# Axes
set xlabel "X [ft]"
set ylabel "Y [ft]" rotate by 0 offset -2,0


# Title (enhanced text, pngcairo-safe)
set title sprintf("%s scheme solving {/Symbol \266}U/{/Symbol \266}t = {/Symbol a} ( {/Symbol \266}^2U/{/Symbol \266}x^2 + {/Symbol \266}^2U/{/Symbol \266}y^2 ) at t = %.2f hr", scheme, tlabel)


# Geometry
set size ratio -1
set tics out
set grid back

# ---- pm3d setup (MATLAB contourf equivalent) ----
set view map
set pm3d map
set pm3d interpolate 0,0   # IMPORTANT: do not smear values

# ---- Color control (FORCED, no autoscale) ----
unset autoscale cb
unset cbrange
set cbrange [0:200]

# Jet-like colormap
set palette rgbformulae 33,13,10

# Colorbar
set colorbox
set cbtics (0, 50, 100, 150, 200)
set cblabel "Temperature" rotate by 0 offset 5,0

unset key

# ---- Plot ----
# Works with XYZ data + blank line between rows
splot datafile using 1:2:3 with pm3d notitle
