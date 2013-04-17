set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417155546625899BkZm9Y/density/density6.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417155546625899BkZm9Y/density/density_custom2.dat" using 1:2 title "Density" with lines lw 3
