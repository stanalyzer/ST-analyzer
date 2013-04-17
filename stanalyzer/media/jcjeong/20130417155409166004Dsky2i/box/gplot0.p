set terminal png
set xlabel 'ps/frame'
set ylabel 'System size'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417155409166004Dsky2i/box/box0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417155409166004Dsky2i/box/box_size.dat" using 1:2 title "X-axis" with linespoints, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417155409166004Dsky2i/box/box_size.dat" using 1:3 title "Y-axis" with linespoints, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130417155409166004Dsky2i/box/box_size.dat" using 1:4 title "Z-axis" with linespoints
