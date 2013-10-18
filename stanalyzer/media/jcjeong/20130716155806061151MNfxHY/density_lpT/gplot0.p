set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716155806061151MNfxHY/density_lpT/density_lpT0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716155806061151MNfxHY/density_lpT/density_lipid_tail.dat" using 1:2 title "Density" with lines lw 3
