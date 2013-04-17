set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417155546625899BkZm9Y/density/density3.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417155546625899BkZm9Y/density/density_water.dat" using 1:2 title "Density" with lines lw 3
