set terminal png
set xlabel 'ps/frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130130114545896684kk2Izt/helix_tilt/helix_tilt0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130130114545896684kk2Izt/helix_tilt/output.dat" using 1:2 title "Helix Tilt" with linespoints
