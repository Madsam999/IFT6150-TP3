/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-B.c                            */
/* Auteur  : Samuel Fournier & Alexandre Toutant        */
/* Date    :                                            */
/* version :                                            */ 
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

#define NAME_IMG_OUT1 "photograph_original_A"
#define NAME_IMG_OUT2 "photograph_degraded_withoutNoise_A" 
#define NAME_IMG_OUT3 "photograph_restored_withoutNoise_A"  
#define NAME_IMG_OUT4 "photograph_degraded_withNoise_A" 
#define NAME_IMG_OUT5 "photograph_restored_withNoise_A"  

void addmatrix(float** resultReal, float** resultComplex, float** mat1Real, float** mat1Complex, float** mat2Real, float** mat2Complex, int width, int height) {
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            resultReal[i][j] = mat1Real[i][j] + mat2Real[i][j];
            resultComplex[i][j] = mat1Complex[i][j] + mat2Complex[i][j];
        }
    }
}

void buildLowPassFilter(float** filterReal, int L, int width, int height) {

    float amplitude = 1 / SQUARE(L);

    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            if((i < L / 2 || i >= height - L / 2) && (j < L / 2 || j >= width - L / 2)) {
                filterReal[i][j] = amplitude;
            }
            else {
                filterReal[i][j] = 0.f;
            }
        }
    }
}

void landweber(float** fReal, float** fComplex, float** gReal, 
               float** gComplex, float** filterReal, float** filterComplex, 
               float** imageReal, int width, int height, int nbIteration,
               float alpha) {
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            fReal[i][j] = gReal[i][j];
        }
    }

    float** intermed1Real = fmatrix_allocate_2d(height, width);
    float** intermed1Complex = fmatrix_allocate_2d(height, width);
    float** intermed2Real = fmatrix_allocate_2d(height, width);
    float** intermed2Complex = fmatrix_allocate_2d(height, width);

    for(int k = 0; k < nbIteration; k++) {
        FFTDD(fReal, fComplex, height, width);

        MultMatrix(intermed1Real, intermed1Complex, filterReal, filterComplex, fReal, fComplex, height, width);

        IFFTDD(fReal, fComplex, height, width);

        IFFTDD(intermed1Real, intermed1Complex, height, width);

        Mult(intermed1Real, -1.0f, height, width);
        Mult(intermed1Complex, -1.0f, height, width);

        addmatrix(intermed2Real, intermed2Complex, intermed1Real, intermed1Complex, gReal, gComplex, width, height);

        FFTDD(intermed2Real, intermed2Complex, height, width);

        MultMatrix(intermed1Real, intermed1Complex, intermed2Real, intermed2Complex, filterReal, filterComplex, height, width);

        IFFTDD(intermed1Real, intermed1Complex, height, width);

        Mult(intermed1Real, alpha, height, width);
        Mult(intermed1Complex, alpha, height, width);

        addmatrix(fReal, fComplex, intermed1Real, intermed1Complex, fReal, fComplex, height, width);       
    }

    SaveImagePgm("temp", fReal, height, width);

    free_fmatrix_2d(intermed1Complex);
    free_fmatrix_2d(intermed1Real);
    free_fmatrix_2d(intermed2Complex);
    free_fmatrix_2d(intermed2Real);
}

int main(int argc, char* argv) {
    int nbIterations;
    int width, height;
    float var, isnr;
    int sizeFilter;

    float** imageReal;
    float** imageComplex;

    imageReal = LoadImagePgm(NAME_IMG_IN, &height, &width);
    imageComplex = fmatrix_allocate_2d(height, width);
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            imageComplex[i][j] = 0.f;
        }
    }

    FFTDD(imageReal, imageComplex, height, width);

    printf("Entrez la largeur du filtre passe bas: ");
    scanf("%d", &sizeFilter);

    float** lowPassFilterReal = fmatrix_allocate_2d(height, width);
    float** lowPassFilterComplex = fmatrix_allocate_2d(height, width);

    float amplitude = 1.0 / SQUARE(sizeFilter);

    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            if((i < sizeFilter / 2.0 || i >= height - sizeFilter / 2.0) && (j < sizeFilter / 2.0 || j >= width - sizeFilter / 2.0)) {
                lowPassFilterReal[i][j] = amplitude;
            }
            else {
                lowPassFilterReal[i][j] = 0.f;
            }
            lowPassFilterComplex[i][j] = 0.f;
        }
    }

    FFTDD(lowPassFilterReal, lowPassFilterComplex, height, width);

    float** gReal = fmatrix_allocate_2d(height, width);
    float** gComplex = fmatrix_allocate_2d(height, width);

    MultMatrix(gReal, gComplex, lowPassFilterReal, lowPassFilterComplex, imageReal, imageComplex, height, width);

    IFFTDD(gReal, gComplex, height, width);
    IFFTDD(imageReal, imageComplex, height, width);
    SaveImagePgm(NAME_IMG_OUT1, imageReal, height, width);
    SaveImagePgm(NAME_IMG_OUT2, gReal, height, width);

    printf("Entrez le nombre d'itérations pour LANDWEBER: ");
    scanf("%d", &nbIterations);

    float** fReal = fmatrix_allocate_2d(height, width);
    float** fComplex = fmatrix_allocate_2d(height, width);

    landweber(fReal, fComplex, gReal, gComplex, lowPassFilterReal, lowPassFilterComplex, imageReal, width, height, nbIterations, 1);

    printf("Entrez la variance du bruit : ");
    scanf("%f",&var);
	
    /* Ajouter du bruit a l'image floue : g = g + bruit = image + flou + bruit (add_gaussian_noise) */
 
    /*******************************************************/
	/* restorer l'image g (elle contient du flou ainsi que */
	/* du bruit) avec LANDWEBER.                           */
	/* N'oubliez pas d'afficher le ISNR a chaque iteration */
	/*******************************************************/
		
    /*Sauvegarde des images */

    /*Liberation memoire*/
  
    /*retour sans probleme*/ 
    printf("\n C'est fini ... \n\n\n");
    return 0; 
}