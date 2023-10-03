#include<time.h>
#include<math.h>
#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"subf.h"
#include"VTSW.h"
#include"data.h"


int main()
{

struct type_donneesc param;
float **x,**y,**xv,**yv,**vol;
float **T0,**T1;
float **U;     // U definie au centre des facettes est-ouest en x,yv
float **V;     // V definie au centre des facettes nord-sud en xv,y
float **Fadv;
float **Fdiff;
float dt,tf;
int l,N;

float **A;   //matrice A
float *B;    //Vecteur second membre
float **AB;  //matrice bande
float *X;    //vecteur qui contiendra la solution de SOR
int LB;      //Largeur de la bande;

printf("simulation started\n");

param=read_datac();
int tailleMat=(param.nx*param.ny);
//cout<<"dans le C : "<<param.nx<<" "<<param.ny<<" "<<param.Lx<<" "<<param.Ly<<" "<<param.U0<<" "<<param.D<<" "<<param.Ti<<" "<<param.Tg<<" "<<param.Tb<<" "<<param.tf<<" "<<param.Nout<<" "<<param.CFL<<" "<<param.R<<endl;
// Allocation dynamique
    x=(float**)malloc((param.nx+1)*sizeof(float*));
    for (int i=0;i<param.nx+1;i++) {x[i]=(float*)malloc((param.ny+1)*sizeof(float));}
    y=(float**)malloc((param.nx+1)*sizeof(float*));
    for (int i=0;i<param.nx+1;i++) {y[i]=(float*)malloc((param.ny+1)*sizeof(float));}
    xv=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {xv[i]=(float*)malloc((param.ny)*sizeof(float));}
    yv=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {yv[i]=(float*)malloc((param.ny)*sizeof(float));}
    vol=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {vol[i]=(float*)malloc((param.ny)*sizeof(float));}
    T0=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {T0[i]=(float*)malloc((param.ny)*sizeof(float));}
    T1=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {T1[i]=(float*)malloc((param.ny)*sizeof(float));}
    U=(float**)malloc((param.nx+1)*sizeof(float*));
    for (int i=0;i<param.nx+1;i++) {U[i]=(float*)malloc((param.ny)*sizeof(float));}
    V=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {V[i]=(float*)malloc((param.ny+1)*sizeof(float));}
    Fadv=(float**)malloc((param.nx)*sizeof(float *));
    for (int i=0;i<param.nx;i++) {Fadv[i]=(float*)malloc((param.ny)*sizeof(float));}
    Fdiff=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {Fdiff[i]=(float*)malloc((param.ny)*sizeof(float));}

meshc(param,x,y,xv,yv,vol);

initial_conditionc(param,xv,yv,x,y,T0,U,V);

VTSWriterc(0.0,0,param.nx+1,param.ny+1,x,y,T0,U,V,"ini");

dt= calc_dtc(param.nx,param.ny,x,y,U,param.CFL,V,param.D,param.R);

N=(int)(param.tf/dt);

printf("N=%d\n",N);

printf("Nout=%d\n",param.Nout);

// Computation with the explicit scheme
if(param.i_solver==0)
  {
  printf("Simulation with the explicit solver...\n");
  printf("time advancing...\n");
  for (l=1;l<N;l++)
    {
    printf("Iteration l=%d...\n",l);
    calc_flux_advc(param,x,y,xv,yv,U,V,T0,Fadv);
    calc_flux_diffc(param,x,y,xv,yv,T0,Fdiff);
    advance_timec(param,dt,vol,Fadv,Fdiff,T0,T1);
    if((l%param.Nout)==0)
      {
      VTSWriterc((float)(l)*dt,l,param.nx+1,param.ny+1,x,y,T1,U,V,"int");
      }
    }
  VTSWriterc((float)(N)*dt,N,param.nx+1,param.ny+1,x,y,T1,U,V,"end");
  }

// Computation with the implicit scheme (Gauss method)
if((param.i_solver)==1)
  {
  printf("Simulation with the implicit solver Gauss...\n");

  A=(float**)malloc((tailleMat)*sizeof(float *));
  for (int i=0;i<tailleMat;i++){A[i]=(float*)malloc((tailleMat)*sizeof(float));}
  if (A==NULL) {printf("Allocation failed");}

  B=(float*)malloc((tailleMat)*sizeof(float*));

  printf("not implemented\n");
   // (to be completed)

  free(A);
  free(B);
  }


// Computation with the implicit scheme (SOR method)
if((param.i_solver)==2)
  {
  printf("Simulation with the implicit solver SOR...\n");

  A=(float**)malloc((tailleMat)*sizeof(float *));
  for (int i=0;i<tailleMat;i++){A[i]=(float*)malloc((tailleMat)*sizeof(float));}
  if (A==NULL) {printf("Allocation failed");}

  B=(float*)malloc((tailleMat)*sizeof(float*));
  LB=(2*param.nx)+1;

  AB=(float**)malloc((tailleMat)*sizeof(float *));
  for (int i=0;i<tailleMat;i++){AB[i]=(float*)malloc((LB)*sizeof(float));}
  if (AB==NULL) {printf("Allocation failed for AB");}
  X=(float*)malloc((tailleMat)*sizeof(float*));

  // (to be completed)
  printf("not implemented\n");


  free(A);
  free(B);
  free(AB);
  free(X);
  }

free(x);
free(y);
free(vol);
free(T0);
free(T1);
free(U);
free(V);
free(Fadv);
free(Fdiff);

printf("Simulation done\n");
printf("dt = %e \n",dt);
return 0;

}


