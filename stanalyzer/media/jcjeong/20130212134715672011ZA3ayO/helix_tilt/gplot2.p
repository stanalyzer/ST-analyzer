set terminal png
set xlabel 'ps/frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130212134715672011ZA3ayO/helix_tilt/helix_tilt2.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130212134715672011ZA3ayO/helix_tilt/output2.dat" using 1:2 title "Helix Tilt" with linespoints
