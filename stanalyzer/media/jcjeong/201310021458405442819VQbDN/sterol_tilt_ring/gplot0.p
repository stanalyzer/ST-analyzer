set terminal png
set encoding iso_8859_1
set xlabel 'Degree'
set ylabel 'Probability'
set title 'Sterol Tilt'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201310021458405442819VQbDN/sterol_tilt_ring/sterol_tilt_ring0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201310021458405442819VQbDN/sterol_tilt_ring/all_dist_sterol_ring_tilt.dat" using 1:2 title "Sterol Ring Tilt" with lines lw 3
