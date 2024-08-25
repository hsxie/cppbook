set term post color solid enh
set output 'planet_orbit.ps'
set multiplot
set origin 0.0,0.0
set size 0.5,1.0
set xlabel 'x'
set ylabel 'y'
plot "data" u 2:3 title 'Euler', "data" u 5:6 title 'Leapfrog'
set origin 0.5,0.0
set xlabel 't'
set ylabel 'total energy'
plot "data" u 1:4 title 'Euler', "data" u 1:7 title 'Leapfrog'