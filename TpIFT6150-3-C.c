/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-C.c                            */
/* Auteur  : Samuel Fournier & Alexandre Toutant        */
/* Date    : 2025-11-05                                 */
/* version : 1.2 (Landweber seul)                       */
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "FonctionDemo3.h"

#define NAME_IMG_IN   "photograph"
#define NAME_IMG_OUT1 "photograph_original_C"
#define NAME_IMG_OUT2 "photograph_degraded_C"
#define NAME_IMG_OUT3 "photograph_restaured_C"


void addmatrix(float** rR, float** rI, float** aR, float** aI, float** bR, float** bI, int W, int H){
    for(int i=0;i<H;i++) for(int j=0;j<W;j++){ rR[i][j]=aR[i][j]+bR[i][j]; rI[i][j]=aI[i][j]+bI[i][j]; }
}

/* PSF boîte LxL centrée (wrap-around) en spatial, amplitude 1/L^2 */
void buildLowPassFilter(float** h, int L, int W, int H){
    for(int i=0;i<H;i++) for(int j=0;j<W;j++) h[i][j]=0.f;
    float amp=1.0f/SQUARE(L); int half=L/2;
    for(int di=0;di<L;di++) for(int dj=0;dj<L;dj++){
        int ii=(di-half+H)%H, jj=(dj-half+W)%W; h[ii][jj]=amp;
    }
}

/* ISNR = 10 log10( sum(f-g)^2 / sum(f-u)^2 ) */
double compute_isnr(float** f, float** g, float** u, int H, int W){
    double num=0.0, den=0.0;
    for(int i=0;i<H;i++) for(int j=0;j<W;j++){ double e1=f[i][j]-g[i][j], e2=f[i][j]-u[i][j]; num+=e1*e1; den+=e2*e2; }
    return 10.0*log10(num/den);
}

float computeDiff() {
    float diff = 0.f;

    return diff;
}

/* Landweber: f^{k+1} = f^k + alpha * H* F(g - H f^k) */
void landweber(float** fR, float** fI, float** gR, float** gI,float** Hre, float** Him, float** fTrue,int W, int H, int iters, float alpha){
    for(int i=0;i<H;i++){
      for(int j=0;j<W;j++){
        fR[i][j]=gR[i][j]; fI[i][j]=0.f; gI[i][j]=0.f;
      }
    }

    float** HcR = Hre;
    float** HcI = fmatrix_allocate_2d(H,W);
    for(int i=0;i<H;i++) for(int j=0;j<W;j++) HcI[i][j] = -Him[i][j];

    float** t1R=fmatrix_allocate_2d(H,W), **t1I=fmatrix_allocate_2d(H,W);
    float** t2R=fmatrix_allocate_2d(H,W), **t2I=fmatrix_allocate_2d(H,W);

    for(int k=0;k<iters;k++){
        FFTDD(fR,fI,H,W);
        MultMatrix(t1R,t1I,Hre,Him,fR,fI,H,W);
        IFFTDD(fR,fI,H,W);
        IFFTDD(t1R,t1I,H,W);

        Mult(t1R,-1.0f,H,W);
        Mult(t1I,-1.0f,H,W);
        addmatrix(t2R,t2I,t1R,t1I,gR,gI,W,H);

        FFTDD(t2R,t2I,H,W);
        MultMatrix(t1R,t1I,HcR,HcI,t2R,t2I,H,W);
        IFFTDD(t1R,t1I,H,W);

        Mult(t1R,alpha,H,W); Mult(t1I,alpha,H,W);
        addmatrix(fR,fI,t1R,t1I,fR,fI,W,H);

        double isnr = compute_isnr(fTrue,gR,fR,H,W);
        printf("Landweber iter %d: ISNR = %.3f dB\n", k+1, isnr);
    }

    free_fmatrix_2d(HcI);
    free_fmatrix_2d(t1I); free_fmatrix_2d(t1R);
    free_fmatrix_2d(t2I); free_fmatrix_2d(t2R);
}

void figueiredoNowak() {
    int k = 0;
    // 1. 25 Iterations of Landweber

    // 2. While not 3x10^-5 do 
    // 2.1 1 Iteration of landweber
    // 2.2 Decomposition
}

int main(int argc, char** argv){
    int H, W, L, iters;
    float var;

    printf("Entrez la taille du filtre passe bas: ");
    scanf("%d", &L);
    printf("Entrez la variance du bruit: ");
    scanf("%f", &var);
    printf("Entrez le nombre d'iterations (Landweber): ");
    scanf("%d", &iters);

    // Lecture de l'image originale
    float** f = LoadImagePgm(NAME_IMG_IN, &H, &W);
    float** fI = fmatrix_allocate_2d(H, W);
    for(int i = 0; i < H; i++) for(int j = 0; j < W; j++) fI[i][j] = 0.f;

    //Construction du filtre passe-bas /
    float** Hre = fmatrix_allocate_2d(H, W);
    float** Him = fmatrix_allocate_2d(H, W);
    buildLowPassFilter(Hre, L, W, H);
    FFTDD(Hre, Him, H, W);  // passage en domaine fréquentiel

    // Simulation du flou
    FFTDD(f, fI, H, W);
    float** gR = fmatrix_allocate_2d(H, W);
    float** gI = fmatrix_allocate_2d(H, W);
    MultMatrix(gR, gI, Hre, Him, f, fI, H, W);
    IFFTDD(gR, gI, H, W);
    IFFTDD(f, fI, H, W);
    SaveImagePgm(NAME_IMG_OUT1, f, H, W);

    //add_gaussian_noise(gR, H, W, var);

    // sauvegarde de l'image dégradée
    float** gView = fmatrix_allocate_2d(H, W);
    copy(gView, gR, H, W);
    Recal(gView, H, W);
    SaveImagePgm(NAME_IMG_OUT2, gView, H, W);
    free_fmatrix_2d(gView);

    // Landweber
    float** uR = fmatrix_allocate_2d(H, W);
    float** uI = fmatrix_allocate_2d(H, W);
    for(int i = 0; i < H; i++) for(int j = 0; j < W; j++) uI[i][j] = 0.f;

    landweber(uR, uI, gR, gI, Hre, Him, f, W, H, iters, 1.0f);

    /*******************************************************************************************
    *
    * Faire 6 itération de l'algo de  Figueiredo et Nowak, devrait donné un ISNR autour de 3.17
    * voir derniere page de l'énoncé
    *
    *******************************************************************************************/

    // Image restaurée
    float** uView = fmatrix_allocate_2d(H, W);
    copy(uView, uR, H, W);
    Recal(uView, H, W);
    SaveImagePgm(NAME_IMG_OUT3, uView, H, W);
    free_fmatrix_2d(uView);

    printf("\nC'est fini ...\n\n");
    return 0;
}