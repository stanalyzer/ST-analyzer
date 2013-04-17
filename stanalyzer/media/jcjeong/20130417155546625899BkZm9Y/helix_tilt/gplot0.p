set terminal png
set xlabel 'ps/frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417155546625899BkZm9Y/helix_tilt/helix_tilt0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417155546625899BkZm9Y/helix_tilt/helix_tilt.dat" using 1:2 title "Helix Tilt" with linespoints
