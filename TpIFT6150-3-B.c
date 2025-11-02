/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-A.c                            */
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

#define NAME_IMG_OUT1 "photograph_bruite_B"
#define NAME_IMG_OUT2 "photograph_debruite_B"


int main(int argc,char** argv){
  int i,j,k,l;
  int length,width;
  int nbLevels;
	
	float threshold;
	float var;   
	
	printf("Entrez la variance du bruit: ");
  scanf("%f",&var);
	
	printf("Entrez le nombre de niveaux a traiter : ");
  scanf("%d",&nbLevels);
		
	printf("Entrez le seuil : ");
  scanf("%f",&threshold);

	/* ouvrir l'image d'entree */

	/* ajouter du bruit a l'image d'entree (add_gaussian_noise) */
	
 	/* debruiter l'image en seuillant les coefficients de Haar */

	/* afficher l'ISBN */
 
 	/* sauvegarder les images */
 
  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
