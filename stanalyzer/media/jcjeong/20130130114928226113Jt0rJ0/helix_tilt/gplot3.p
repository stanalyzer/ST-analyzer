set terminal png
set xlabel 'ps/frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130130114928226113Jt0rJ0/helix_tilt/helix_tilt3.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130130114928226113Jt0rJ0/helix_tilt/output3.dat" using 1:2 title "Helix Tilt" with linespoints
