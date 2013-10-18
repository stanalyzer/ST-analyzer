set terminal png
set encoding iso_8859_1
set xlabel 'Degree'
set ylabel 'Probability'
set title 'Sterol Tilt'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130718141250144697wrSElg/sterol_tilt_tail/sterol_tilt_tail0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130718141250144697wrSElg/sterol_tilt_tail/all_dist_sterol_tail_tilt.dat" using 1:2 title "Sterol Tail Tilt" with lines lw 3
