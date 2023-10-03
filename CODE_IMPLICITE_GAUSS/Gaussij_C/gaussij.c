#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<math.h>

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
