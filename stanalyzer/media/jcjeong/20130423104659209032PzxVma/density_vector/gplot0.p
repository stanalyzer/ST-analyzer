set terminal png
set xlabel 'degree'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130423104659209032PzxVma/density_vector/density_vector0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130423104659209032PzxVma/density_vector/density_water_dipole.dat" using 1:2 title "Density" with lines lw 3
