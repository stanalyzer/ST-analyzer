set terminal png
set xlabel 'ps/frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130212135022242593tVJDbu/helix_tilt/helix_tilt20.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130212135022242593tVJDbu/helix_tilt/output20.dat" using 1:2 title "Helix Tilt" with linespoints
