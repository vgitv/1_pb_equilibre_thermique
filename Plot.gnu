# réinitialisation des paramètres
reset

# fichier de sortie
set term postscript eps size 3.5,2.62 enhanced color font 'Helvetica,12'
set output "graphes/u.eps"
set encoding utf8

# paramètres
set title "Équation de transport à vitesse constante"
set grid
set xlabel "x"
set ylabel "y"
#set size ratio -1 # repère orthonormé
#set xrange[-2:2]
#set yrange[-1.5:1.5]

# tracé
plot "sorties/approx.dat" u 1:2     lc rgb "#008000" lw 1 title "approximation",\
     "sorties/sol.dat"    u 1:2 w l lc rgb "#FF4500" lw 1 title "solution"

# affichage écran
set term wxt enhanced
replot
