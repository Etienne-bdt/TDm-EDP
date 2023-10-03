#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//*******************************************************************************************************************
void  VTSWriterc(float Time, int Step,int nx,int ny,float **x,float **y,float **T,float **U,float **V,char opt[3])  {
  //----------------------------------------------------------------------------
  //  Time    : Reel, temps physique
  //Step      : Entier, pas de temps = numero dans le nom de fichier,           #
  //          si Step<0 on ecrit !                                             #
  //          dans sol_exacte.vts, entier                                      #
  // nx      : Entier, nombre des noeuds en direction x                         #
  // ny      : Entier, nombre des noeuds en direction y                         #
  // x       : Tableau reel (de taille nx,ny) des abscisses des noeuds          #
  //            des  volumes                                                     #
  // y       : Tableau reel (de taille nx,ny) des ordonnees des noeuds          #
  //           des  volumes                                                     #
  // T       : Tableaux reel (de taille nx-1 par ny-1) des valeurs a tracer     #
  //          (valeurs au centre des volumes de controle)                      #
  // U       : Tableaux reel (de taille nx   par ny-1) des valeurs a tracer     #
  //           (valeurs au centre des facettes de normale + ou -x               #
  // V       : Tableaux reel (de taille nx-1 par ny  ) des valeurs a tracer     #
  //           (valeurs au centre des facettes de normale + ou -y               #
  // opt     : Variable de type chaine des characteres qui doit prendre         #
  //         l'une des valeurs suivantes :                                    #
  //           - 'ini' pour le premier appel a VTSWriter                      #
  //           - 'int' pour un appel standard a VTSWriter                     #
  //           - 'end' pour le dernier appel a VTSWriter                      #
  //-----------------------------------------------------------------------------#
  FILE * fp;
  FILE * vp;
  char num2char[100];
  char FileName[200];
  int i,j;

  //  --- Ecriture d'un fichier temporel au format paraview  ---

  // nom du fichier : sol_Step.vts, où Step est le numéro de l'itération en temps
  snprintf(FileName, 100, "sol_%d.vts", Step);
  // ouverture du fichier FileName
    fp = fopen (FileName,"w");
  // ecriture du header
    fprintf( fp, "<?xml version=\"1.0\"?>\n" );
    fprintf( fp, "<VTKFile type=\"StructuredGrid\">\n" );
    fprintf( fp, "<StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0,nx-1,0,ny-1,0,0 );
    fprintf( fp, "<Piece Extent=\"%d %d %d %d %d %d\">\n", 0,nx-1,0,ny-1,0,0 );
  // ecriture des coordonnees du maillage

       fprintf( fp, "<Points>\n" );
       fprintf( fp, "<DataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n" );

       for (j = 0; j < ny; j++){
          for (i = 0; i < nx; i++){
             fprintf( fp, "%f %f %f ", x[i][j], y[i][j], 0. );
          }
          fprintf( fp, "\n" );
       }


       fprintf( fp, "</Points>\n" );
  // Début d'écriture des scalaires

       fprintf( fp, "<CellData Scalars=\"Temperature, U, V\">\n" );

  // Ecriture de la température

       fprintf( fp, "<DataArray type=\"Float32\" Name=\"Temp, K\"/>\n" );

       for (j = 0; j < ny-1; j++){
          for (i = 0; i < nx-1; i++){
             fprintf( fp, "%f ", T[i][j] );
          }
          fprintf( fp, "\n" );
       }

   // écriture de la vitesse u

       fprintf( fp, "<DataArray type=\"Float32\" Name=\"Vitesse u, m/s\"/>\n" );

       for (j = 0; j < ny-1; j++){
          for (i = 0; i < nx-1; i++){
             fprintf( fp, "%f ", (U[i+1][j] + U[i][j])/2.0 );
          }
          fprintf( fp, "\n" );
       }

    // Ecriture de la vitesse v

       fprintf( fp, "<DataArray type=\"Float32\" Name=\"Vitesse v, m/s\"/>\n" );

       for (j = 0; j < ny-1; j++){
          for (i = 0; i < nx-1; i++){
             fprintf( fp, "%f ", (V[i][j+1] + V[i][j])/2.0 );
          }
          fprintf( fp, "\n" );
       }

   // Fin d'écriture des scalaires
       fprintf( fp, "</CellData>\n" );

   // Ecriture du footer
       fprintf( fp, "</Piece>\n" );
       fprintf( fp, "</StructuredGrid>\n" );
       fprintf( fp, "</VTKFile>\n" );

   // Fermeture du fichier FileName
       fclose (fp);

   // Remplissage du fichier "Collection" determinant l evolution temporelle -

   if (strcmp( opt, "ini") == 0 )
   { vp= fopen ("sol.pvd","w");
       fprintf( fp, "<?xml version=\"1.0\"?>\n" );
       fprintf( fp, "<VTKFile type=\"Collection\" version=\"0.1\" format=\"ascii\">\n" );
       fprintf( fp, "<Collection>\n" ); }
   else{
       vp= fopen ("sol.pvd","a");
   }
   if (Step >= 0)
       {fprintf( vp, "<DataSet timestep=\"%f\" group=\"\" part=\"\" file=\"%s\" />\n", Time, FileName ) ;}
   if (strcmp( opt, "end") == 0 )
       {fprintf( fp, "</Collection>\n" );
        fprintf( fp, "</VTKFile>\n" );
   }
                 fclose(vp);

  }

