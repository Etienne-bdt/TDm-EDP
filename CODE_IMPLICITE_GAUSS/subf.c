#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include"data.h"



//*********************************
struct type_donneesc read_datac(){

   struct    type_donneesc param;
   FILE *fp;
   char str[1000], str_copy[1000];
   int iter;
   int read_int;
   float read_float;

   printf("reading data...\n");
   fp = fopen("data.txt", "r");
   // lit les donnees à partir du fichier
      iter = 0;

       while ( fgets (str, 1000, fp) != NULL ) {

           strncpy(str_copy, str, 10);

           if ( iter == 0 || iter == 1 || iter == 10 || iter == 13 || iter == 14 || iter == 15  ){
               read_int = atoi(str_copy);
           }
           else{
               read_float = atof(str_copy);
           }

           switch (iter){

               case 0:
                   param.nx = read_int;
                   break;

               case 1:
                   param.ny = read_int;
                   break;

               case 2:
                   param.Lx = read_float;
                   break;

               case 3:
                   param.Ly = read_float;
                   break;

               case 4:
                   param.U0 = read_float;
                   break;

               case 5:
                   param.D = read_float;
                   break;

               case 6:
                   param.Ti = read_float;
                   break;

               case 7:
                   param.Tg = read_float;
                   break;

               case 8:
                   param.Tb = read_float;
                   break;

               case 9:
                   param.tf = read_float;
                   break;

               case 10:
                   param.Nout = read_int;
                   break;

               case 11:
                   param.CFL = read_float;
                   break;

               case 12:
                   param.R = read_float;
                   break;

               case 13:
                   param.i_mesh = read_int;
                   break;

               case 14:
                   param.i_vit = read_int;
                   break;

               case 15:
                   param.i_solver = read_int;
                   break;

               case 16:
                   param.W = read_float;
                   break;

               case 17:
                   param.R0 = read_float;
                   break;

           }

           iter++;

       }

       fclose(fp);

       return param;

   }


//********************************************************************************************
void meshc(struct type_donneesc param,float **x,float **y,float **xv,float **yv,float **vol) {

int i,j;
printf("creating mesh...\n");

// Uniform mesh (streching-compressing flow + rotating flow)
//for (i=0;i<param.nx+1;i++){
//    for (j=0;j<param.ny+1;j++){
//     x[i][j]=-param.Lx+2.*param.Lx*(float)(i)/(float)(param.nx);
//     y[i][j]=-param.Ly+2.*param.Ly*(float)(j)/(float)(param.ny);
//    }
//}

// Uniform mesh (parabolic flow)
for (i=0;i<param.nx+1;i++){
    for (j=0;j<param.ny+1;j++){
     x[i][j]=param.Lx*(float)(i)/(float)(param.nx);
     y[i][j]=param.Ly*(float)(j)/(float)(param.ny);
    }
}


for (i=0;i<param.nx;i++){
    for( j=0;j<param.ny;j++){
          xv[i][j]=0.5*(x[i][j]+x[i+1][j]);
          yv[i][j]=0.5*(y[i][j]+y[i][j+1]);
          vol[i][j]=(x[i+1][j]-x[i][j])*(y[i][j+1]-y[i][j]);
    }
}

if(param.i_mesh == 1){
    for(i=0;i<param.nx+1;i++){
        for(j=0;j<param.ny+1;j++){
            x[i][j] = (pow(x[i][j],2))/param.Lx;
            y[i][j] = param.Ly*(1-cos(M_PI*y[i][j]/(2*param.Ly)));
        }
    }
    for(i=0;i<param.nx;i++){
        for(j=0;j<param.ny;j++){
            xv[i][j]=0.5*(x[i][j]+x[i+1][j]);
            yv[i][j]=0.5*(y[i][j]+y[i][j+1]);
            vol[i][j]=(x[i+1][j]-x[i][j])*(y[i][j+1]-y[i][j]);
        }
    }
};

}




//****************************************************************************************************************************
void initial_conditionc(struct type_donneesc param,float **xv,float **yv,float **x,float **y,float **T0,float **U,float **V){

    int i,j;

    printf("initial condition...\n");

    for (i=0;i<param.nx;i++){
        for (j=0;j<param.ny;j++){
            T0[i][j]=param.Ti;
        }
    }

    if( param.i_vit == 0){
    for (i=0;i<param.nx+1;i++){
        for (j=0;j<param.ny;j++){
            U[i][j] = param.U0 ;
        }
    }
    }
    else{
        for(i=0;i<param.nx+1;i++){
            for(j=0;j<param.ny;j++){
                U[i][j]= 6*param.U0*(y[i][j]/(2*param.Ly))*(1-(y[i][j]/(2*param.Ly)));
            }
        }
    };
    for (i=0;i<param.nx;i++){
        for (j=0;j<param.ny+1;j++){
            V[i][j] = 0.;
        }
    }


      }


//*****************************************************************************************************
float calc_dtc(int nx, int ny,float **x,float **y,float **U,float CFL,float **V,float D,float R) {

 float dx,dy;
 float dt,dt_loc;
 int  i,j;

 printf("computing dt...\n");

 dt = 1.e8;

 for( j=0; j<ny; j++){
   for (i=0; i<nx; i++){
     dx = x[i+1][j]-x[i][j];
     dy = y[i][j+1]-y[i][j];
     dt_loc = 1.0/(fabs(U[i][j])/(CFL*dx)+fabs(V[i][j])/(CFL*dy)+D*(1.0/(dx*dx)+1.0/(dy*dy))/R);
    if (dt_loc < dt) { dt = dt_loc; }
   }
 }

 printf("dt=%e\n",dt);

return dt;

   }



//************************************************************************************************************************************
void calc_flux_advc(struct type_donneesc param,float **x,float **y,float **xv,float **yv,float **U,float **V,float **T0,float **Fadv)
{

    float **Fadv_o, **Fadv_e, **Fadv_n, **Fadv_s;
    float T_amont;
    int i,j,k;
    Fadv_e=(float**)malloc((param.nx)*sizeof(float *));
    for (int i=0;i<param.nx;i++) {Fadv_e[i]=(float*)malloc((param.ny)*sizeof(float));}
    Fadv_o=(float**)malloc((param.nx)*sizeof(float *));
    for (int i=0;i<param.nx;i++) {Fadv_o[i]=(float*)malloc((param.ny)*sizeof(float));}
    Fadv_n=(float**)malloc((param.nx)*sizeof(float *));
    for (int i=0;i<param.nx;i++) {Fadv_n[i]=(float*)malloc((param.ny)*sizeof(float));}
    Fadv_s=(float**)malloc((param.nx)*sizeof(float *));
    for (int i=0;i<param.nx;i++) {Fadv_s[i]=(float*)malloc((param.ny)*sizeof(float));}

    // Advection term

   // East flux
   // In the domain
    for (i=0;i<param.nx-1;i++){
        for (j=0;j<param.ny;j++){
            if (U[i+1][j]>=0.0){
                T_amont = T0[i][j];}
            else{
                T_amont = T0[i+1][j];
            }
            Fadv_e[i][j] = -1.0*U[i+1][j]*T_amont*(y[i+1][j+1]-y[i+1][j]);
        }
    }
   // BC
    i=param.nx-1;
    for (j=0;j<param.ny;j++){
        if (U[i+1][j]>=0.0){
            T_amont = T0[i][j];}
        else{
            printf("BC not implemented 01 \n");
            break;
        }
        Fadv_e[i][j] = -1.0*U[i+1][j]*T_amont*(y[i+1][j+1]-y[i+1][j]);
    }
   //West flux
   //In the domain
    for (i=1;i<param.nx;i++){
        for (j=0;j<param.ny;j++){
            Fadv_o[i][j] = -1.0*Fadv_e[i-1][j];
        }
    }
   // BC
    i=0;
    for (j=0;j<param.ny;j++){
        if (U[i][j]<=0.0){
            T_amont = T0[i][j];}
        else{
            T_amont = param.Tg;
        }
        Fadv_o[i][j] = U[i][j]*T_amont*(y[i][j+1]-y[i][j]);
    }
   // North flux
   // In the domain
    for (i=0;i<param.nx;i++){
        for (j=0;j<param.ny-1;j++){
            if (V[i][j+1]>=0.0){
                T_amont = T0[i][j];}
            else{
                T_amont = T0[i][j+1];
            }
            Fadv_n[i][j] = -1.0*V[i][j+1]*T_amont*(x[i+1][j+1]-x[i][j+1]);
        }
    }
   // BC
    j=param.ny-1;
    for (i=0;i<param.nx;i++){
        if (V[i][j+1]>=0.0){
            T_amont = T0[i][j];}
       else{
            T_amont = 0.0;
       }
       Fadv_n[i][j] = -1.0*V[i][j+1]*T_amont*(x[i+1][j+1]-x[i][j+1]);
    }
   //South flux
   //In the domain
    for (i=0;i<param.nx;i++){
        for (j=1;j<param.ny;j++){
            Fadv_s[i][j] = -1.0*Fadv_n[i][j-1];
        }
    }
   // BC
    j=0;
    for (i=0;i<param.nx;i++){
        if (V[i][j]<=0.){
            T_amont = T0[i][j];}
        else{
            T_amont = param.Tb;
        }
        Fadv_s[i][j] = V[i][j]*T_amont*(x[i+1][j]-x[i][j]);

    }

   //c) Total summation

    for (i=0;i<param.nx;i++){
        for (j=0;j<param.ny;j++){
            Fadv[i][j] = Fadv_o[i][j] + Fadv_e[i][j] + Fadv_s[i][j] + Fadv_n[i][j];
        }
    }

   }

//*********************************************************************************************************************************
void calc_flux_diffc(struct type_donneesc param,float **x, float **y,float **xv,float **yv,float **T0, float **Fdiff) {

float  **Fdiff_o, **Fdiff_e, **Fdiff_s, **Fdiff_n;
int  i,j,k;
Fdiff_e=(float**)malloc((param.nx)*sizeof(float *));
for (int i=0;i<param.nx;i++) {Fdiff_e[i]=(float*)malloc((param.ny)*sizeof(float));}
Fdiff_o=(float**)malloc((param.nx)*sizeof(float *));
for (int i=0;i<param.nx;i++) {Fdiff_o[i]=(float*)malloc((param.ny)*sizeof(float));}
Fdiff_n=(float**)malloc((param.nx)*sizeof(float *));
for (int i=0;i<param.nx;i++) {Fdiff_n[i]=(float*)malloc((param.ny)*sizeof(float));}
Fdiff_s=(float**)malloc((param.nx)*sizeof(float *));
for (int i=0;i<param.nx;i++) {Fdiff_s[i]=(float*)malloc((param.ny)*sizeof(float));}

// Diffusive term
// West flux
// In the domain
for (i=1;i<param.nx;i++){
    for (j=0;j<param.ny;j++){
        Fdiff_o[i][j] = -1.0*param.D*(T0[i][j]-T0[i-1][j])/(xv[i][j]-xv[i-1][j])*(y[i][j+1]-y[i][j]);
    }
}
//BC
i=0;
for (j=0;j<param.ny;j++){
    Fdiff_o[i][j] = -1.0*param.D*(T0[i][j]-param.Tg)/(2.*(xv[i][j]-x[i][j]))*(y[i][j+1]-y[i][j]);
}
//East flux
// In the domain
for (i=0;i<param.nx-1;i++){
    for (j=0;j<param.ny;j++){
        Fdiff_e[i][j] = -1.0*Fdiff_o[i+1][j];
    }
}
//BC
i=param.nx-1;
for (j=0;j<param.ny;j++){
    Fdiff_e[i][j] = -1.0*Fdiff_o[i][j];
}
// South flux
//In the domain
for (i=0;i<param.nx;i++){
    for (j=1;j<param.ny;j++){
        Fdiff_s[i][j] = -1.0*param.D*(T0[i][j]-T0[i][j-1])/(yv[i][j]-yv[i][j-1])*(x[i+1][j]-x[i][j]);
    }
}
//BC
j=0;
for (i=0;i<param.nx;i++){
    Fdiff_s[i][j] = -1.0*param.D*(T0[i][j]-param.Tb)/(2.*(yv[i][j]-y[i][j]))*(x[i+1][j]-x[i][j]);
}
// North flux
// In the domain
for (i=0;i<param.nx;i++){
    for (j=0;j<param.ny-1;j++){
        Fdiff_n[i][j] = -1.0*Fdiff_s[i][j+1];
    }
}
//BC
j=param.ny-1;
for (i=0;i<param.nx;i++){
    Fdiff_n[i][j] = 0.0;
}

//c) Total summation

for (i=0;i<param.nx;i++){
    for (j=0;j<param.ny;j++){
        Fdiff[i][j] = Fdiff_o[i][j] + Fdiff_e[i][j] + Fdiff_s[i][j] + Fdiff_n[i][j];
    }
}

}

//*********************************************************************************************************************
void  advance_timec(struct type_donneesc param,float dt, float **vol,float **Fadv,float **Fdiff,float **T0,float **T1){

int i,j;

for (i=0;i<param.nx;i++)
{
    for (j=0;j<param.ny;j++)
    {
        T1[i][j] = T0[i][j] + dt/vol[i][j]*(Fadv[i][j] + Fdiff[i][j]);
    }
}
for (i=0;i<param.nx;i++)
{
    for  (j=0;j<param.ny;j++)
    {
        T0[i][j]=T1[i][j];
    }
}
}

//*********************************************
void creation_A(struct type_donneesc param,int NA, float dt, float **x,float **y,float **xv,float **yv,float **vol,float **A)
{

int i,j,k;
float a,b,c,d,e;
float dx,dy;

for(i=0;i<NA;i++){
    for(j=0;j<NA;j++){
        A[i][j]=0;
    }
}

for(j=0;j<param.ny;j++){
    for(i=0;i<param.nx;i++){
        dx= x[i+1][j]-x[i][j];
        dy= y[i][j+1]-y[i][j];

        //Dans ce code, on se propose de créer la matrice A en commencant par ses coefficients. 
        //Il est à noter que les coefficients a,b,c,d,e sont reliés par une relation linéaire de sorte que
        //a = - (b+c+d+e) dans le cas général. Dans certains cas cette relation devient fausse si l'un des termes
        //est pris en compte dans B, pour pallier à ce problème, nous allons le calculer de la même façon que dans B, 
        //le prendre en compte dans a puis l'annuler.

        switch (i)
        {
        case 0:
            d = -(dt*param.D/vol[i][j])*(dy)/(xv[i][j]*2);
            e = -(dt*param.D/vol[i][j])*(dy)/(x[i+1][j]/2);
            break;
        
        case 1:
            d = -(dt*param.D/vol[i][j])*(dy)/(xv[i][j]-xv[i-1][j]);
            e = -(dt*param.D/vol[i][j])*(dy)/(xv[i-1][j]*2);
            break;

        default:
            d = -(dt*param.D/vol[i][j])*(dy)/(xv[i][j]-xv[i-1][j]);
            e = -(dt*param.D/vol[i][j])*(dy)/(xv[i-1][j]-xv[i-2][j]);
            break;
        }

        switch (j)
        {
        case 0:
            b = -(dt*param.D/vol[i][j])*(dx)/(yv[i][j]*2);
            c = -(dt*param.D/vol[i][j])*(dx)/(y[i][j+1]/2);
            break;
        
        case 1:
            b = -(dt*param.D/vol[i][j])*(dx)/(yv[i][j]-yv[i][j-1]);
            c = -(dt*param.D/vol[i][j])*(dx)/(yv[i][j-1]*2);
            break;

        default:
            b = -(dt*param.D/vol[i][j])*(dx)/(yv[i][j]-yv[i][j-1]);
            c = -(dt*param.D/vol[i][j])*(dx)/(yv[i][j-1]-yv[i][j-2]);
            break;
        }

        if(i==0&&j==0){
            //Coin bas gauche
            a = -(b+c+d+e);
            c = 0;//Pris en compte dans B
            e = 0; //Pris en compte dans B
        }
        else if(i==0&&j!=0&&j!=param.ny-1){//Frontière gauche
            a = -(b+c+d+e);
            e = 0;//Pris en compte dans B
        }
        else if(i==0&&j==param.ny-1){//Coin Haut Gauche
            b =0;//Annulé
            a =  -(b+c+d+e);
            e = 0;//Pris en compte dans B
        }
        else if(i!=0&&j==param.ny-1&&i!=param.nx-1){
            //Frontière Haute
            b = 0;//Annulé
            a = -(b+c+d+e);
        }
        else if(i==param.nx-1&&j==param.ny-1){
            //Coin Haut Droit
            b = 0;//Annulé
            d = 0;//Annulé
            e = 0;//Annulé
            a = -(b+c+d+e);
        
        }
        else if(j!=0&&i==param.nx-1&&j!=param.ny-1){
            //Frontière Droite
            d = 0;//Annulé
            e = 0;//Annulé
            a = -(b+c+d+e);
        }
        else if(i==param.ny-1&&j==0){
            //Coin Bas Droit
            d = 0;//Annulé
            e = 0;//Annulé
            a = -(b+c+d+e);
            c = 0;//Pris en compte dans B
        }
        else if(i!=0&&j==0){
            //Frontière basse
            a = -(b+c+d+e);
            c = 0; //Pris en compte dans B 
        }
        else{
            //Cas général
            a = -(b+c+d+e);
        }
        
        //Assignation des valeurs de A
        k = i+j*param.nx;
        
        if (k>0){
            A[k][k-1] = e;
        }
        if (k<param.nx*param.ny){
            A[k][k+1] = d;
        }
        if(k>param.nx-1){ 
            A[k][k-param.nx] = c;
        }
        if(k<param.ny*param.nx-1){
            A[k][k+param.nx] = b;
        }
        A[k][k]=1+a;
        }
    }
}


//**********************************


void creation_B(struct type_donneesc param, int NA, float dt, float **x, float **y,float **xv,float **yv,float **vol, float **Fadv, float **T0, float *B)
{
int i,j,k;
float dx,dy;
float a,b,c,d,e;
for (j=0;j<param.ny;j++){
    for(i=0;i<param.nx;i++){
        dx= x[i+1][j]-x[i][j];
        dy= y[i][j+1]-y[i][j];
        k = i+j*param.nx;
        if(i==0&&j!=param.ny-1&&j!=0){
            //Frontiere Gauche
            e = -(dt*param.D/vol[i][j])*(dy)/((x[i+1][j]/2)); 
            B[k] = T0[i][j]+((dt/vol[i][j])*Fadv[i][j]-e*param.Tg);
        }

        else if(i!=0 && j==0&&i!=param.nx-1){
            //Frontiere Bas
            c = -(dt*param.D/vol[i][j])*(dx)/(y[i][j+1]/2);
            B[k] = T0[i][j]+((dt/vol[i][j])*Fadv[i][j]-c*param.Tb);
        }
        else if (i==0 && j==param.ny-1){
            //Coin haut-gauche
            e = -(dt*param.D/vol[i][j])*(dy)/((x[i+1][j]/2)); 
            B[k] = T0[i][j]+((dt/vol[i][j])*Fadv[i][j]-e*param.Tg);
        }
        else if (i==0 && j==0){
            //Coin bas-gauche
            c = -(dt*param.D/vol[i][j])*(dx)/(y[i][j+1]/2);
            e = -(dt*param.D/vol[i][j])*(dy)/((x[i+1][j]/2)); 
            B[k] = T0[i][j]+((dt/vol[i][j])*Fadv[i][j]-c*param.Tb - e*param.Tg);
        }
        else if (i==param.nx-1&&j==0){
            //Coin bas-droite
            c = -(dt*param.D/vol[i][j])*(dx)/(y[i][j+1]/2);
            B[k] = T0[i][j]+((dt/vol[i][j])*Fadv[i][j]-c*param.Tb);
        }
        else {
            B[k] = T0[i][j]+((dt/vol[i][j])*Fadv[i][j]);
            
        }
        }
    }
}



//*****************************************

void miseajour_T(struct type_donneesc param,float **T0,float **T1,float *B)
{
int i,j,k;

for (j=0;j<param.ny;j++){
    for(i=0;i<param.nx;i++){
        k = i+j*param.nx;
        T0[i][j] = B[k];
        T1[i][j] = B[k];
    }
}

}

//***********************************

//*******************************
//
//             RESOLUTION D'UN SYSTEME LINEAIRE
//
// METHODE : Methode de Gauss.
//
// LANGAGE : C
//
//  MODE D'UTILISATION :  GAUSSIJ (LV,A,B)
//
//  Donnees : A  Coefficients de la matrice, variable a deux dimensions
//              dont les valeurs numeriques doivent etre fournies
//              conformement au schema ci-dessous.
//           B  Termes du 2eme membre, variable a un indice.
//               A la sortie de GAUSS, la solution se trouve dans B.
//            LV Ordre de la matrice A.
//
//        |                                                       |
//        | A(1,1)       *                  *                    *|
//        |  *                                                    |
//        |                                                       |
//        |                                                       |
//        |                                                       |
//        |                                                       |
//        |                                                       |
//        |                                                       |
//        |                                                       |
//        |                                                       |
//        |                                                       |
//        |  *                                    *       A(LV,LV)|
//
void GAUSSIJ(int LV, float **A, float *B){
    
    float **ATEMP;
    int I,J,K;
    float DIVB,BI,AKI;
    ATEMP=(float**)malloc((LV)*sizeof(float *));
    for (int I=0;I<LV;I++) {ATEMP[I]=(float*)malloc((LV)*sizeof(float));}

    for(I=0; I<LV;I++){
        for(J=0;J<LV;J++){
            ATEMP[I][J]=A[I][J];
        }
    }
    for(I=0; I<LV;I++){
        DIVB=1.0/ATEMP[I][I];
        B[I]=B[I]*DIVB;
        BI=B[I];
        for(J=LV-1; J>I-1;J--){
            ATEMP[I][J]=ATEMP[I][J]*DIVB;
        }
        for(K=0;K<I;K++){
            AKI=ATEMP[K][I];
            B[K]=B[K]-AKI*BI;
            for(J=LV-1;J>I-1;J--){
                ATEMP[K][J]=ATEMP[K][J]-AKI*ATEMP[I][J];
            }
        }
        for(K=I+1;K<LV;K++){
            AKI=ATEMP[K][I];
            B[K]=B[K]-AKI*BI;
            for(J=LV-1;J>I-1;J--){
                ATEMP[K][J]=ATEMP[K][J]-AKI*ATEMP[I][J];
            }
        }
    }

    return;

}
//

//***********************************
//***********************************
//
//           PASSAGE D'UNE MATRICE PLEINE A UNE MATRICE BANDE
//
// LANGAGE : C
//
//  Donnees : Cette routine cree une matrice bande AB de dimension (NA,LB)
//            a partir d'une matrice pleine A de dimension (NA,NA) telle que:
//           NA  Ordre de la matrice A
//           LB Largeur de la bande = 2 MX + 1;
//            MX est appelee "demi largeur de bande"
//
//             La matrice A est de la forme:
//
//             | A(1,1)       *                  *                    *|
//             |  *                                                    |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |  *                                    *       A(NA,NA)|
//
//
//             Les coefficients de la matrice AB correspondent a:
//
//AB(1,1) ... | AB(1,MX+1)   ...    AB(1,LB)                            |
//             |                                                         |
//             | AB(2,MX)      .....    AB(2,LB)                         |
//             |                                                         |
//             |           *       .....       *                         |
//             |                                                         |
//             |               *       .....       *                     |
//             |                                                         |
//             |                   AB(I,1) ... AB(I,MX+1) ... AB(I,LB)   |
//             |                                                         |
//             |                         *       .....       *           |
//             |                                                         |
//             |                             *       .....       *       |
//             |                                                         |
//             |                               AB(NA,1)  ...  AB(NA,MX+1)| ...  AB(NA,LB)

void create_A_band(int MX,int NA,float **A, float **AB){

    int i,j;
//printf("hala esaie1\n");
    for(i=0;i<NA;i++){
        for(j=0;j<2*MX+1;j++){
            AB[i][j] = 0.0;
        }
    }
//printf("hala esaie2\n");
    for(i=0;i<MX+1;i++){
//printf("hala esaieplz\n");
        for(j=0;j<NA-(i);j++){
//printf("hala esaie10\n");
            AB[j][(MX)+(i)] = A[j][j+(i)];
//printf("hala esaie11\n");
            AB[NA-1-j][(MX)-(i)] = A[NA-1-j][NA-1-j-(i)];
//printf("hala esaie12\n");
        }
   //printf("hala esaie3\n"); 
}
//printf("hala esaie4\n");
       return;

       }
//***************************

//************************************
//************************************
//           RESOLUTION D'UN SYSTEME LINEAIRE
//
//  METHODE : Methode de Sur-Relaxation Successive SOR.
//
//  LANGAGE : C
//
//  Donnees : A  Coefficients de la matrice bande, variable a deux dimensions
//              dont les valeurs numeriques doivent etre fournies conformement
//               au schema ci-dessous.
//              (Les points exterieurs a la bande ne sont pas pris en compte
//              lors du calcul).
//            B  Termes du 2eme membre, variable a un indice
//            X  A la sortie de SOR, la solution se trouve dans X
//           N  Ordre de la matrice A
//            LB Largeur de la bande = 2 MX + 1
//            MX est appelee "demi largeur de bande"
//            R0 precision souhaitee
//           W  parametre de sur-relaxation (0<W<2)
//
//            |                                                       |
// A(1,1) ... | A(1,MX+1)   ...    A(1,LB)                            |
//            |                                                       |
//            | A(2,MX)      .....    A(2,LB)                         |
//            |                                                       |
//            |           *       .....       *                       |
//            |                                                       |
//            |               *       .....       *                   |
//            |                                                       |
//            |                   A(I,1) ... A(I,MX+1) ... A(I,LB)    |
//            |                                                       |
//            |                         *       .....       *         |
//            |                                                       |
//            |                             *       .....       *     |
//            |                                                       |
//            |                                 A(N,1)  ...  A(N,MX+1)| ...  A(N,LB)

void SOR(int MX,int N,float **A,float *B,float R0,float W, float *X){

    int  i,cc,dd,compteur;
    float *Y = malloc(sizeof(int)*(2*MX+N));
    float *Y0= malloc(sizeof(int)*(2*MX+N));
    float ECART,C,D,Z,Z0;
    

    for(i=0;i<2*MX+N;i++){
        Y[i]=0.0;
        Y0[i]=0.0;
    }
    ECART=1.0;
    compteur = 0;
    do  {
        compteur = compteur + 1;
        for(i=MX;i<MX+N;i++){
          // -------> calcul de la somme inférieure appelée C
            C=0.0;
            for(cc=0;cc<MX;cc++){
                C=C+Y[i-1-cc]*A[i-MX][MX-cc-1];
            }
          // -------> calcul de la somme supérieure appelée D
            D=0.0;
            for(dd=0;dd<MX;dd++){
                D=D+Y0[i+dd+1]*A[i-MX][MX+1+dd];
            }
          // -------> calcul de X,k+1
            Y[i]=(1.-W)*Y0[i]+(W/A[i-MX][MX])*(B[i-MX]-C-D);
        }

        // -------> critère d'arrêt
        Z0=0.0;
        Z=0.0;
        for(i=MX;i<MX+N;i++){
            Z=Z+Y[i]*Y[i];
            Z0=Z0+Y0[i]*Y0[i];
        }
      
        ECART=fabs(pow(Z,0.5)-pow(Z0,0.5));
        
        for(i=MX; i<MX+N;i++){
            Y0[i]=Y[i];
        }

    }
    while ((ECART>R0)&&(compteur<=10000));

    if (compteur>=10000) {
        printf("Probleme de convergence du solveur SOR\n");
    }

      // -------> écriture dans X

    for(i=MX;i<MX+N;i++){
        X[i-MX]=Y[i];
    }

    return;

}
//********************

