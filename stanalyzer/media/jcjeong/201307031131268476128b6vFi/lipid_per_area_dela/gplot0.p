set terminal png
set xlabel 'Time (ps)'
set ylabel 'Area per Lipid (A^2)'
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201307031131268476128b6vFi/lipid_per_area_dela/lipid_per_area_dela0.png'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/201307031131268476128b6vFi/lipid_per_area_dela/delaunay_lipid_per_area.dat" using 1:2 title "custom thicknesss" with lines lw 3
