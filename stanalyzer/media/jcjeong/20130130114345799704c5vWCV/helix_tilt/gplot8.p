set terminal png
set xlabel 'ps/frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130130114345799704c5vWCV/helix_tilt/helix_tilt8.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130130114345799704c5vWCV/helix_tilt/output8.dat" using 1:2 title "Helix Tilt" with linespoints
