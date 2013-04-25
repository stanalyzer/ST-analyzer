set terminal png
set xlabel 'range'
set ylabel '<cos0>'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130425094911790068d2Xsk3/vector_selection/vector_selection0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130425094911790068d2Xsk3/vector_selection/density_vector.dat" using 1:2 title "Density" with lines lw 3
