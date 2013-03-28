set terminal png
set xlabel 'range'
set ylabel 'density'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130328132603481668w7rhmO/density/density0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130328132603481668w7rhmO/density/output.dat" using 1:2 title "Density" with lines lw 3
