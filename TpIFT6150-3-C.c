/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-D.c                            */
/* Auteur  : Francois Destrempes                        */
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

#define NAME_IMG_OUT1 "photograph_original_C"
#define NAME_IMG_OUT2 "photograph_degraded_C" 
#define NAME_IMG_OUT3 "photograph_restaured_C"  

int main(int argc,char** argv)
 {
  int nb_iterations;
  int i,j,k,l;
  int length,width;
  float var;
  int size_filter;
  		
  printf("Entrez la taille du filtre passe bas : ");
  scanf("%d",&size_filter);	

  printf("Entrez la variance du bruit : ");
  scanf("%f",&var);
 
  printf("Entrez le nombre d'itérations (LANDWEBER): ");
  scanf("%d",&nb_iterations);

  /* ouvrir l'image d'entree */

  /* ajouter du flou et du bruit a l'image d'entree (add_gaussian_noise) */
 	
	/* 1er etape : deconvolution avec 'nb_iterations' de LANDWEBER */
	
	/* 2e etape : Filtrage dans le domaine des ondelettes */
    
  /* Sauvegarde des matrices sous forme d'image pgm */

  /* Liberation memoire pour les matrices */

  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
