set terminal png
set xlabel 'Frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130111120741798733BKSQVH/helix_tilt/helix_tilt.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130111120741798733BKSQVH/helix_tilt/output.dat" using 1:2 title "Helix Tilt" with linespoints
