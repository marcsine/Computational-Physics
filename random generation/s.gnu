set term png
set title 'Normalized Histogram of ρ(x)'
set grid

set output 'fig1.png'
plot 'his.dat' i 0 u 1:2:3 w yerrorbars t'Histogram',(12/(pi*(2*pi**2-3)))*x**2*(sin(x))**2 t"ρ(x)"

set title 'Normalized Histogram – Exponential Distribution'
set output 'fig2.png'
plot 'his.dat' i 1 u 1:2:3 w yerrorbars t'Histogram',pi*exp(-pi*x) t"λe^{−λx}"