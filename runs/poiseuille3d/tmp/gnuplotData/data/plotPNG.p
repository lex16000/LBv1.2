set key right
set terminal png
set output './tmp/gnuplotData/centerVelocity.png'
set xlabel ''
set ylabel ''
plot './tmp/gnuplotData/data/centerVelocity.dat' u 1:2 w l t 'analytical', './tmp/gnuplotData/data/centerVelocity.dat' u 1:3 w l t 'numerical'
