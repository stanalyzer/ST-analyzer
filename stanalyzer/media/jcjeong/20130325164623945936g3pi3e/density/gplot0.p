set terminal png
set xlabel 'ps/frame'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130325164623945936g3pi3e/density/density0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130325164623945936g3pi3e/density/output.dat" using 2:3 title "Density" with linespoints
