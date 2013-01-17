set terminal png
set xlabel 'Frame'
set ylabel 'System size'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130114094836629779ECkPZO/box/box.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130114094836629779ECkPZO/box/output.dat" using 1:2 title "X-axis" with linespoints, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130114094836629779ECkPZO/box/output.dat" using 1:3 title "Y-axis" with linespoints, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130114094836629779ECkPZO/box/output.dat" using 1:4 title "Z-axis" with linespoints
