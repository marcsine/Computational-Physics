set term png
set grid

set output "figP1.png"
set yrange [-0.5:6]

plot 'res11.dat' i 0 u 2:3 w l t"#1",\
'res11.dat' i 1 u 2:3 w l t"#2",\
'res11.dat' i 2 u 2:3 w l t"#3",\
'res11.dat' i 3 u 2:3 w l t"#4"


#esta que viene ahora no la piden

set output "other.png"
set yrange [0:4.5]

plot 'res22.dat' i 0 u 1:2 w l notitle