Nom du programme : prog_exp.exe

Date de création : 12/09/23

Auteurs : OUEDRAOGO - MARRO - Thomas Bonometti


Objectif : Programme (en C) de resolution d'une equation d'advection-diffusion 2D. Methode des volumes finis. Schema Euler explicite en temps, amont pour l'advection, centre pour la diffusion.

Fichier d'entrée: "data.txt"

Fichiers résultats : "sol_000000.vts", "sol_000001.vts", etc,

Pour compiler le programme :
Dans un terminal, se placer dans le dossier contenant les sources du programme (*.c, Makefile), puis taper: 
make

Pour exécuter le programme : 
Dans un terminal, se placer dans le dossier contenant le fichier executable (prog_exp.exe) et le fichier d'entree (data.txt), puis taper:
./prog_exp.exe 

Pour visualiser les résultats:
Dans un terminal, se placer dans le dossier contenant les fichiers de resultat (*.vts), puis taper:
paraview &
