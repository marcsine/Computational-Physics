set term png
set output 'f.png'
set yrange [-1.5:2.5]
set grid

set term pngcairo font ",10"
set xtics font ",8"

set tmargin 5

set multiplot layout 1, 2 #font ",14"
set label 1 "Eigenfunctions (x in L, Energy in ħ/(m L^{2})" at screen 0.5,0.98 center font ",14" offset 0,-1

set title '50 grid points'

plot 'res.dat' i 0 u 1:2 w l t '1', \
'res.dat' i 1 u 1:2 w l t '2', \
'res.dat' i 2 u 1:2 w l t '3', \
'res.dat' i 3 u 1:2 w l t '4'

set title '640 grid points'

plot 'res.dat' i 4 u 1:2 w l t '1', \
'res.dat' i 5 u 1:2 w l t '2', \
'res.dat' i 6 u 1:2 w l t '3', \
'res.dat' i 7 u 1:2 w l t '4'

unset multiplot
