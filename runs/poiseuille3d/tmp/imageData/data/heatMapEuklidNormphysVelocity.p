if (strstrt(GPVAL_TERMINALS, 'jpeg') > 0) {set terminal jpeg size 1920,1080font ",25"
set output './tmp/imageData/_EuklidNorm(physVelocity)iT0013112.jpeg'
} else {set terminal png size 1920,1080font ",25"
set output './tmp/imageData/_EuklidNorm(physVelocity)iT0013112.png'
}set pm3d map
unset key
set size ratio -1
set size 0.925,1.0
set xtics out
set ytics out
set xtics nomirror
set ytics nomirror
set pm3d interpolate 0,0
set colorbox vertical user origin 0.85,0.1 size 0.025 ,0.8
set xlabel "x-axis in m "
set ylabel "y-axis in m "
set cblabel offset 0.5 "EuklidNorm(physVelocity) in m/s"
set autoscale fix
set palette defined ( 0 "blue", 1 "green", 2 "yellow", 3 "orange", 4 "red" )
splot './tmp/imageData/data/EuklidNorm(physVelocity).matrix' u ($1*0.0338889+-0.0555556):($2*0.0338889+-2.22045e-16):3 matrix with pm3d
