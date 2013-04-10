set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201304091444083596076IJInM/density/density0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201304091444083596076IJInM/density/density_lipid_head.dat" using 1:2 title "Density" with lines lw 3
