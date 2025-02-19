set term png

set output 'fig1.png'

set xlabel 'E'
set xrange[0:5]
set yrange[-10:10]

set grid

plot 'p3-1920-res.dat' i 0 u 1:2 w l t 'F(E)', '' i 0 u 1:3 w l t 'dF(E)/dE'
unset xrange
unset yrange

######################################################

set output 'fig2.png'
set title 'Newton-Raphson for different initial values'
set xlabel 'Iterations'

plot 'conv_NR.dat' i 0 u 1:2 w l t '0.2','' i 1 u 1:2 w l t '1.5','' i 2 u 1:2 w l t '0.3'

######################################################

set output 'fig3.png'
set title 'Approximating F(E)'
set yrange [-10000:70000]
set xlabel 'E'
plot 'P3-1920-res3_n34.dat' u 1:2 w p t'34 Grid Points', 'P3-1920-res3_n420.dat' u 1:3 w p t'420 Grid Points', \
'P3-1920-res3_n420.dat' u 1:4 w l t'Exact'
