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

void saveImage(char* outFile, float** image, int height, int width) {
    float** outImage = fmatrix_allocate_2d(height, width);
    copy(outImage, image, height, width);
    Recal(outImage, height, width);
    SaveImagePgm(outFile, outImage, height, width);
    free_fmatrix_2d(outImage);
}

void addmatrix(float** rR, float** rI, float** aR, float** aI, float** bR, float** bI, int W, int H){
    for(int i=0;i<H;i++) for(int j=0;j<W;j++){ rR[i][j]=aR[i][j]+bR[i][j]; rI[i][j]=aI[i][j]+bI[i][j]; }
}

float diffSquared(float** mat1, float** mat2, int width, int height) {
    float sum = 0.f;
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            sum += SQUARE(mat1[i][j] - mat2[i][j]);
        }
    }
    return sum;
}

float computeISNR(float** fReal, float** gReal, float** imageReal, int width, int height) {
    float isnr;
    isnr = 10 * log10f(diffSquared(imageReal, gReal, width, height) / diffSquared(imageReal, fReal, width, height));
    return isnr;
}

void landweber(float** fReal, float** fComplex, float** gReal, 
               float** gComplex, float** filterReal, float** filterComplex, 
               float** imageReal, int width, int height, int nbIteration,
               float alpha) {

    float** intermed1Real = fmatrix_allocate_2d(height, width);
    float** intermed1Complex = fmatrix_allocate_2d(height, width);
    float** intermed2Real = fmatrix_allocate_2d(height, width);
    float** intermed2Complex = fmatrix_allocate_2d(height, width);

    for(int k = 0; k < nbIteration; k++) {
        double isnr = computeISNR(fReal,gReal,imageReal,width,height);
        printf("Landweber iter %d: ISNR = %f\n", k+1, isnr);
        FFTDD(fReal, fComplex, height, width);
        FFTDD(gReal, gComplex, height, width);

        MultMatrix(intermed1Real, intermed1Complex, filterReal, filterComplex, fReal, fComplex, height, width);

        IFFTDD(fReal, fComplex, height, width);
        IFFTDD(gReal, gComplex, height, width);

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

    free_fmatrix_2d(intermed1Complex);
    free_fmatrix_2d(intermed1Real);
    free_fmatrix_2d(intermed2Complex);
    free_fmatrix_2d(intermed2Real);
}

int main(int argc, char** argv){
    int maxIters;
    int width, height;
    float var, alpha;
    int sizeFilter;
    int nbLevels = 3;
    int nbIters;

    float** imageReal;
    float** imageComplex;
    float** zeros;

    // Load initial image and initialize imaginary matrix
    imageReal = LoadImagePgm(NAME_IMG_IN, &height, &width);
    imageComplex = fmatrix_allocate_2d(height, width);
    // Helper matrix
    zeros = fmatrix_allocate_2d(height, width);
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            imageComplex[i][j] = 0.f;
            zeros[i][j] = 0.f;
        }
    }

    printf("Entrez la largeur du filtre passe bas: ");
    scanf("%d", &sizeFilter);
    printf("Entrez la variance du bruit (sigma^2): ");
    scanf("%f",&var);
    printf("Entrez le nombre d'itérations initiales pour Landweber: ");
    scanf("%d", &nbIters);

    // Creation of sizeFilter x sizeFilter filter
    float** hReal = fmatrix_allocate_2d(height, width);
    float** hComplex = fmatrix_allocate_2d(height, width);

    float amplitude = 1.0 / SQUARE(sizeFilter);

    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            if((i < sizeFilter / 2.0 || i >= height - sizeFilter / 2.0) && (j < sizeFilter / 2.0 || j >= width - sizeFilter / 2.0)) {
                hReal[i][j] = amplitude;
            }
            else {
                hReal[i][j] = 0.f;
            }
            hComplex[i][j] = 0.f;
        }
    }

    // Fourier Transform
    FFTDD(hReal, hComplex, height, width);
    FFTDD(imageReal, imageComplex, height, width);

    // Create result image of f x h (image with blur)
    float** gReal = fmatrix_allocate_2d(height, width);
    float** gComplex = fmatrix_allocate_2d(height, width);

    // Convolution --FFT-> Product
    MultMatrix(gReal, gComplex, imageReal, imageComplex, hReal, hComplex, height, width);

    // Revert Fourier Transform
    IFFTDD(gReal, gComplex, height, width);
    IFFTDD(imageReal, imageComplex, height, width);

    // Add gaussian noise to blurred image: g = f*h + n -> (blurred image with gaussian noise)
    add_gaussian_noise(gReal, height, width, var);

    saveImage(NAME_IMG_OUT1, imageReal, height, width);
    saveImage(NAME_IMG_OUT2, gReal, height, width);

    float** fReal = fmatrix_allocate_2d(height, width);
    float** fComplex = fmatrix_allocate_2d(height, width);

    // Initial call of landweber
    // f_0 = g;
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            fReal[i][j] = gReal[i][j];
        }
    }

    landweber(fReal, fComplex, gReal, gComplex, hReal, hComplex, imageReal, width, height, nbIters, 1.f);

    // Figueiredo & Nowak Algo
    float** fPrev = fmatrix_allocate_2d(height, width);

    char loop = 1;
    float isnr;
    int k = 0;

    // In prep for Haar
    float** haar = fmatrix_allocate_2d(height, width);
    int avgSize = height / powf(2, nbLevels);

    while(loop) {   
        for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                fPrev[i][j] = fReal[i][j];
            }
        }

        // One iteration of landweber
        landweber(fReal, fComplex, gReal, gComplex, hReal, hComplex, imageReal, width, height, 1, 1.f);

        haar2D_complete(fReal, haar, nbLevels, height, width);

        for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                if(i < avgSize && j < avgSize) {
                    continue;
                }
                if(haar[i][j] != 0) {
                    haar[i][j] = fmax(0, SQUARE(haar[i][j]) - 3 * var) / haar[i][j];
                }
            }
        }

        ihaar2D_complete(haar, fReal, nbLevels, height, width);

        // Calculer le ISNR        
        isnr = computeISNR(fReal, gReal, imageReal, width, height);

        printf("%02d - ISNR : %lf\n", k+1, isnr);

        float num = diffSquared(fReal, fPrev, width, height);
        float denum = diffSquared(fPrev, zeros, width, height);

        loop = (num / denum) > 0.0001;

        k++;
    }

    Recal(fReal, height, width);
    SaveImagePgm(NAME_IMG_OUT3, fReal, height, width);

    // Retour sans problème
    printf("\nC'est fini ...\n\n");
    return 0;
}