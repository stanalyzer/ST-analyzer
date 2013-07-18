set terminal png
set output '/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/lipid_per_area_voro0.png'
set multiplot layout 2, 1 title 'Area per lipid'
set tmargin 2
set title 'Top Membrane'
set xlabel 'Time (ps)'
set ylabel 'Area per Lipid [{Å}^2]'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/top_ave_voronoi_lipid_per_area1401_1403.dat" using 1:2 title "CHL1" with lines lw 3, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/top_ave_voronoi_lipid_per_area1401_1403.dat" using 1:3 title "POPC" with lines lw 3, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/top_ave_voronoi_lipid_per_area1401_1403.dat" using 1:4 title "POPI" with lines lw 3, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/top_ave_voronoi_lipid_per_area1401_1403.dat" using 1:5 title "POPI13" with lines lw 3, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/top_ave_voronoi_lipid_per_area1401_1403.dat" using 1:6 title "DOPC" with lines lw 3
set title 'Bottom Membrane'
set xlabel 'Time (ps)'
set ylabel 'Area per Lipid [{^305}^2]'
plot "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/btm_ave_voronoi_lipid_per_area1401_1403.dat" using 1:2 title "CHL1" with lines lw 3, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/btm_ave_voronoi_lipid_per_area1401_1403.dat" using 1:3 title "POPC" with lines lw 3, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/btm_ave_voronoi_lipid_per_area1401_1403.dat" using 1:4 title "POPI" with lines lw 3, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/btm_ave_voronoi_lipid_per_area1401_1403.dat" using 1:5 title "POPI13" with lines lw 3, "/home2/jcjeong/project/stanalyzer0/stanalyzer/media/jcjeong/20130716133121871146zdDzT0/lipid_per_area_voro/btm_ave_voronoi_lipid_per_area1401_1403.dat" using 1:6 title "DOPC" with lines lw 3
unset multiplot
set output
