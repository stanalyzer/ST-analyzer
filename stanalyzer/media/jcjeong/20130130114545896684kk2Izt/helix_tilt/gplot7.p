set terminal png
set xlabel 'ps/frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130130114545896684kk2Izt/helix_tilt/helix_tilt7.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130130114545896684kk2Izt/helix_tilt/output7.dat" using 1:2 title "Helix Tilt" with linespoints
