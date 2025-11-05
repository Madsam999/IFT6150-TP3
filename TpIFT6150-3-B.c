/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-A.c                            */
/* Auteur  : Samuel Fournier & Alexandre Toutant        */
/* Date    : 2025-11-05                                 */
/* version : 1.0                                        */
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo3.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/  
/*------------------------------------------------*/
#define NAME_IMG_IN  "photograph"

#define NAME_IMG_OUT1 "photograph_bruite_B"
#define NAME_IMG_OUT2 "photograph_debruite_B"


// Seuillage avec renormalisation
static void threshold_block(float **haar, int r0, int r1, int c0, int c1, float T) {
    const float MAXV = 255.0f;
    const float scale = MAXV / (MAXV - T);
    for (int i = r0; i < r1; i++) {
        for (int j = c0; j < c1; j++) {
            float v = haar[i][j];
            float av = fabsf(v);
            if (av < T) {
                haar[i][j] = 0.0f;
            } else {
                haar[i][j] = copysignf((av - T) * scale, v);
            }
        }
    }
}

// ISNR = 10 log10(sum(f-g)^2 / sum(f-u)^2)
static double compute_isnr(float **f, float **g, float **u, int H, int W) {
    double num = 0.0,den = 0.0;
    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++) {
            double e_noisy  = f[i][j] - g[i][j];
            double e_denois = f[i][j] - u[i][j];
            num+= e_noisy * e_noisy;
            den += e_denois * e_denois;
        }
    return 10.0 * log10(num / den);
}

int main(int argc, char **argv) {
    int length, width, nbLevels;
    float threshold, var;

    printf("Entrez la variance du bruit: ");
    scanf("%f", &var);
    printf("Entrez le nombre de niveaux a traiter: ");
    scanf("%d", &nbLevels);
    printf("Entrez le seuil: ");
    scanf("%f", &threshold);

    float **f = LoadImagePgm(NAME_IMG_IN, &length, &width);
    float **g = fmatrix_allocate_2d(length, width);
    float **haar = fmatrix_allocate_2d(length, width);
    float **u = fmatrix_allocate_2d(length, width);

    // Ajout du bruit
    copy(g, f, length, width);
    add_gaussian_noise(g, length, width, var);

    // Sauvegarde visuelle de g sans l’altérer
    float **gView = fmatrix_allocate_2d(length, width);
    copy(gView, g, length, width);
    Recal(gView, length, width);
    SaveImagePgm(NAME_IMG_OUT1, gView, length, width);
    free_fmatrix_2d(gView);

    // Débruitage par seuillage Haar
    haar2D_complete(g, haar, nbLevels, length, width);

    for (int lev = 1; lev <= nbLevels; lev++) {
        int a = 1 << (lev - 1);
        int rHalf = length / (2* a);
        int rFull = length / a;
        int cHalf = width / (2 * a);
        int cFull = width / a;

        threshold_block(haar, 0, rHalf, cHalf, cFull, threshold);      // HL
        threshold_block(haar, rHalf, rFull, 0, cHalf,  threshold);      // LH
        threshold_block(haar, rHalf, rFull, cHalf,cFull, threshold);  // HH
    }

    ihaar2D_complete(haar, u, nbLevels, length, width);

    double isnr = compute_isnr(f, g, u, length, width);
    printf("ISNR (dB) = %.3f\n", isnr);

    Recal(u, length, width);
    SaveImagePgm(NAME_IMG_OUT2, u, length, width);

    printf("\nC'est fini ...\n\n");
    return 0;
}