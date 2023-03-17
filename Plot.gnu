# réinitialisation des paramètres
reset

# fichier de sortie
set term postscript 0 eps size 3.5,2.62 enhanced color font 'Helvetica,12'
set output "graphes/poisson.eps"
set encoding utf8

# paramètres
set title "Équation de Poisson"
set grid
set xlabel "x"
set ylabel "y"
#set size ratio -1 # repère orthonormé
#set xrange[-2:2]
#set yrange[-1.5:1.5]

# tracé
plot "sorties/approx.dat" u 1:2     lc rgb "#008000" lw 1 title "approximation",\
     "sorties/sol.dat"    u 1:2 w l lc rgb "#FF4500" lw 1 title "solution"



set term postscript 1 eps size 3.5,2.62 enhanced color font 'Helvetica,12'
set output "graphes/eq_thermique.eps"
set encoding utf8

# paramètres
set title "Équation à l'équilibre thermique"
set grid
set xlabel "x"
set ylabel "y"
#set size ratio -1 # repère orthonormé
#set xrange[0:3]
set yrange[0:3]

# tracé
 plot "sorties/approx_n.dat"   u 1:2     lc rgb "#008000" lw 1 title "approx n",\
     "sorties/approx_p.dat"   u 1:2     lc rgb "#FF4500" lw 1 title "approx p",\
     "sorties/approx_psi.dat" u 1:2     lc rgb "#1E90FF" lw 1 title "approx psi"

# affichage écran
set term wxt enhanced
replot
