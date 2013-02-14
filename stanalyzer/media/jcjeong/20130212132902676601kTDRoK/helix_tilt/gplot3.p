set terminal png
set xlabel 'ps/frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130212132902676601kTDRoK/helix_tilt/helix_tilt3.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130212132902676601kTDRoK/helix_tilt/output3.dat" using 1:2 title "Helix Tilt" with linespoints
