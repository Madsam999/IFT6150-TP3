/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-B.c                            */
/* Auteur  :                                            */
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

int main(int argc,char** argv)
 {
  int nb_iterations;
  int length,width;
  float var, isnr;
  int size_filter; /* taille du filtre servant a ajouter du flou a l'image d'entree */

  float** image; /* image d'entree */
  float** g;  	 /* image degradee */
  float** f; 		 /* image restoree */
   
	printf("Entrez la largeur du filtre passe bas : ");
  scanf("%d",&size_filter);

  printf("\nEntrez le nombre d'itérations pour LANDWEBER: ");
  scanf("%d",&nb_iterations);
		
  /* lire l'image d'entree */	
	
  /* Ajouter du flou a l'image d'entree : g = image + flou */



  /*******************************************************/
	/* restorer l'image g (elle contient du flou mais pas  */
	/* de bruit) avec LANDWEBER.                           */
	/* N'oubliez pas d'afficher le ISNR a chaque iteration */
	/*******************************************************/
  
  /*Sauvegarde des images */



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
