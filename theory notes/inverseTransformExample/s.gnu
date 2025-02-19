set term png
set output 'fig.png'

plot 'res.dat' u 1:2:3 w yerrorbars,\
1./2.*sin(x) w l
