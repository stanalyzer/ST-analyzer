set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130422164905442379IS2d6z/density/density0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130422164905442379IS2d6z/density/density_water.dat" using 1:2 title "Density" with lines lw 3
