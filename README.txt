Réalisation de tests et simulations

Activation des tests / simulations
Aller dans main.cpp
Modifier les booléens de lancement des différents tests et simulations pour choisir lesquels seront réalisés
Dans le cas des simulations, spécifier le chemin du maillage à utiliser (data/square.mesh, ou data/square_fine.mesh, pour les simulations autres que le mug)
Dans le cas des tests, il est possible de modifier des paramètres à partir de src/tests.h

Pour lancer une simulation différente de celles pré-implémentées : 
- aller dans ./src/tests.h
- aller dans le dernier test (test_solve_poisson)
- modifier le maillage en entrée, et les diverses fonctions utilisées à l'intérieur, à l'aide de la définition de la fonction (dans src/simu.h)

Compilation et réalisation
Dans le terminal, se placer dans fem2A_base_MailleOkada
- make
- pour lancer les simulations précédemment choisies : ./build/fem2A -s
- pour lancer les test précédemment choisis : ./build/fem2A -t
