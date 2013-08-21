set terminal png
set encoding iso_8859_1
set xlabel 'ps/frame'
set ylabel 'System size'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130820153859392310u5v64G/box/box0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130820153859392310u5v64G/box/box_size.dat" using 1:2 title "X-axis" with lines lw 3, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130820153859392310u5v64G/box/box_size.dat" using 1:3 title "Y-axis" with lines lw 3, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130820153859392310u5v64G/box/box_size.dat" using 1:4 title "Z-axis" with lines lw 3
