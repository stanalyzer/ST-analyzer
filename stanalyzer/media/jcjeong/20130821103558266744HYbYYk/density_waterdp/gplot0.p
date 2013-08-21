set terminal png
set xlabel 'range'
set ylabel '<cos0>'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130821103558266744HYbYYk/density_waterdp/density_waterdp0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130821103558266744HYbYYk/density_waterdp/density_water_dipole.dat" using 1:2 title "Density" with lines lw 3
