set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20131004164356257254aKNqfT/density_water/density_water0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20131004164356257254aKNqfT/density_water/density_water.dat" using 1:2 title "Density" with lines lw 3
