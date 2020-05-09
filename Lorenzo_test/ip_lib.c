/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"


/**
 * Modifiche Lorenzo
 **/

/**
 * TODO: trova il minimo di una matrice ip_mat del canale k passato
 **/
float min(ip_mat * a,unsigned int k){
    
    int i,j;
    int altezza_matrice, larghezza_matrice;
    
    int min =   a->data[k][0][0];
    
    
    altezza_matrice =   a->h;
    larghezza_matrice   =  a->w;

    for(i=0;i<altezza_matrice;i++){
        
        for(j=0;j<larghezza_matrice;j++){
            
            if(min>a->data[k][i][j])
                min=a->data[k][i][j];
        }
    }

    return min;
    
}


/**
  *TODO: trova il massimo della matrice ip_mat del canale k passato 
  **/
float max(ip_mat * a,unsigned int k){

    int i,j;
    int altezza_matrice, larghezza_matrice;
    
    int max=a->data[k][0][0];
    
    
    altezza_matrice =   a->h;
    larghezza_matrice   =  a->w;

    for(i=0;i<altezza_matrice;i++){
        
        for(j=0;j<larghezza_matrice;j++){
            
            if( max  <   a->data[k][i][j] )
                max=a->data[k][i][j];
        }
    }

    return max;
    
}


float mean(ip_mat * a,unsigned int k){

    int i,j;
    int altezza_matrice, larghezza_matrice;
    int count=0;   
    int mean=0;
    
    
    altezza_matrice =   a->h;
    larghezza_matrice   =  a->w;

    for(i=0;i<altezza_matrice;i++){
        for(j=0;j<larghezza_matrice;j++){

            mean+=a->data[k][altezza_matrice][larghezza_matrice];
            count++;
        }
    }

    return (mean/count);
}

/* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */
ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){
   int i,j,c;

   ip_mat *nuova;

    /**
     * asseggno alla struttura nuova i valori dei campi
     * dell'altezza e larghezzza
     * */
    /**creazione delle k matrici + inizializzazione a v**/

   nuova->data=(float ***) malloc(sizeof(unsigned int)*k);


    for(c=0;c<k;c++){

        nuova->data[c]= (float**) malloc(sizeof(float)*h);
    }


    for(c=0;c<k;c++){

        for(i=0;i<h;i++){

            nuova->data[k][i]=(float *) malloc(sizeof(float)*w);
        }
    }


    /*inzializzazione dei valori della matrice*/
    
    for(c=0;c<k;c++){
        for(i=0;i<h;i++){
            for(j=0;j<w;j++){
                nuova-> data[k][i][j]    =   v;
            }
        }
    }

    nuova->h = h;
    nuova->k = k;
    nuova->w = w;

    
    /**
     * riempio i valri di stats per ogni canale;
     **/

    for(c=0;c<k;c++){
        nuova->stat->max=max(nuova,c);
        nuova->stat->min=min(nuova,c);
        nuova->stat->mean=mean(nuova,c);
    
    }



  return nuova;    
}



void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    return cos(2*PI*y2)*sqrt(-2.*log(y1));

}