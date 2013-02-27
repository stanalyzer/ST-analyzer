set terminal png
set xlabel 'ps/frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130215093652750891EyevkZ/helix_tilt/helix_tilt0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130215093652750891EyevkZ/helix_tilt/output.dat" using 1:2 title "Helix Tilt" with linespoints
