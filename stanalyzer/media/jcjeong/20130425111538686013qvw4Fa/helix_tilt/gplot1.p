set terminal png
set xlabel 'ps/frame'
set ylabel 'Tilt Angle'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130425111538686013qvw4Fa/helix_tilt/helix_tilt1.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130425111538686013qvw4Fa/helix_tilt/helix2.dat" using 1:2 title "Helix Tilt" with lines lw 3
