Le code s'organise en 5 fichiers: 

-constantes.f90: contient la précision des réels, et éventuellement d'autres constantes utiles

-functions.f90: contient les conditions initiales, conditions de bord et solutions exactes

-Legendre.f90: contient ce qui est relatif aux polynômes de Legendre: poids/abcisses des quadratures, polynomes de Legendre jusqu'au degré 9,
                norme des polynômes, et calcul de la solution approchée à partir de coordonnées dans la base de Legendre 

-matrix.f90: contient ce qui est relatif à la création des matrices/vecteurs pour le schéma: quadrature des matrices L, M, et N, création des matrices,
                quadrature pour la projection des conditions initiales et de bord pour créer les vecteurs alpha/beta. 

-main.f90: contient l'éxécution du code de manière générale 

Il y a un makefile:
- make pour compiler
- ./run pour éxécuter

Dans gnuplot:
-load 'plot.txt' pour tracer les solutions exacte et approchée au temps final choisi.
-load 'plit_erreur.txt' pour tracer la courbe d'erreur
            