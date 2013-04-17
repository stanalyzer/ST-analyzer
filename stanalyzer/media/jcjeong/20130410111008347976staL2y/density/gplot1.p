set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130410111008347976staL2y/density/density1.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130410111008347976staL2y/density/density_lipid_head.dat" using 1:2 title "Density" with lines lw 3
