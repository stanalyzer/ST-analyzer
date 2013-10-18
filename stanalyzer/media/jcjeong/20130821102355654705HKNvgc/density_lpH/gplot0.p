set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130821102355654705HKNvgc/density_lpH/density_lpH0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130821102355654705HKNvgc/density_lpH/density_lipid_head.dat" using 1:2 title "Density" with lines lw 3
