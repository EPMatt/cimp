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
 * Dato un puntatore a ip_mat e un canale ritorna il valore più piccolo della matrice presente in quel canale
 *TODO: se vuota restituiamo 0 o FLT_MIN?
 **/
float min(ip_mat * a,unsigned int k){
    
    unsigned int i,j;
    unsigned int altezza_matrice, larghezza_matrice;
    /*minimo valore psitivo disponibile */
    float min=FLT_MIN;
    
    if(a){

        min =   a->data[k][0][0];
        
        
        altezza_matrice =   a->h;
        larghezza_matrice   =  a->w;

        for(i=0;i<altezza_matrice;i++){
            
            for(j=0;j<larghezza_matrice;j++){
                
                if(min>a->data[k][i][j])
                    min=a->data[k][i][j];
            }
        }
    }

    return min;
    
}

/**
 * Dato un puntatore a ip_mat e un canale ritorna il valore più grande della matrice presente in quel canale
 *TODO: se vuota restituiamo 0 o FLT_max?
 **/

float max(ip_mat * a,unsigned int k){

    unsigned int i,j;
    unsigned int altezza_matrice, larghezza_matrice;

    float max=FLT_MAX;

    if(a){
        max=a->data[k][0][0];
        
        altezza_matrice =   a->h;
        larghezza_matrice   =  a->w;

        for(i=0;i<altezza_matrice;i++){
            
            for(j=0;j<larghezza_matrice;j++){
                
                if( max  <   a->data[k][i][j] )
                    max=a->data[k][i][j];
            }
        }
    }

    return max;
    
}


/**
 * Dato un puntatore a ip_mat e un canale ritorna il valore medio di quel canale in base al numero di elementi
 * 
 * TODO: Controllare: è giusto fare la media in quel modo?
 **/ 


float mean(ip_mat * a,unsigned int k){

    unsigned int i,j;
    unsigned int altezza_matrice, larghezza_matrice;

    int count=0;   
    float mean=0;
    
    if(a){
        altezza_matrice =   a->h;
        larghezza_matrice   =  a->w;

        for(i=0;i<altezza_matrice;i++){
            for(j=0;j<larghezza_matrice;j++){

                mean+=a->data[k][altezza_matrice][larghezza_matrice];
                count++;
            }
        }
    }

    return (mean/count);
}



/**
 * Calcola stat : data una matrice ip_mat calcola e modifica le statistiche di quel canale
 * 
 **/ 

void calcola_stat(ip_mat *a){
    int l=0;
    
    for(l=0;l<a->k;l++){

        a->stat[l].max=max(a,l);
        a->stat[l].min=min(a,l);
        a->stat[l].mean=mean(a,l);
    }

}



/**
 * Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){
   unsigned int i,j,c;

   ip_mat *nuova;
   
   nuova = (ip_mat*) malloc(sizeof(ip_mat)); 
    /**
     * asseggno alla struttura nuova i valori dei campi
     * dell'altezza e larghezzza
     *
     **/
    nuova->h = h;
    nuova->k = k;
    nuova->w = w;


    /**creazione delle k matrici + inizializzazione a v**/
   
    if(nuova){
        /**
         * Creo la matrice a 3 dimensioni effettiva e i 3 vettori con
         * le statistiche per canale
         **/
        nuova->data=(float ***) malloc(sizeof(float**)*k);
        
        nuova->stat=(stats*) malloc(sizeof(stats)*k);

            /**
             * canale
             **/
            for(c=0;c<k;c++){

                nuova->data[c]= (float**) malloc(sizeof(float*)*h);
            }


            for(c=0;c<k;c++){

                for(i=0;i<h;i++){

                    nuova->data[c][i]=(float *) malloc(sizeof(float)*w);
                }
            }


            /*inzializzazione dei valori della matrice*/
            
            for(c=0;c<k;c++){
                for(i=0;i<h;i++){
                    for(j=0;j<w;j++){
                        nuova-> data[c][i][j]    =   v;
                    }
                }
            }

            
            /**
             * riempio i valri di stats per ogni canale;
             * 
             **/

            for(c=0;c<k;c++){
                nuova->stat[c].max=v;
                nuova->stat[c].min=v;
                nuova->stat[c].mean=mean(nuova,c);
            }
    }

  return nuova;
  

}



float random(float mean,float var){
    return 1.0;
}

/**
 * TODO: random(mean, var); 
 * Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e varianza var 
 * 
 * */

void ip_mat_init_random(ip_mat * t, float mean, float var){
    int i,j,l;

    /**
     * 
     * riempio la matrice del canale i-esimo con valori casuali
     *
     * */

    for(l=0;l<t->k;l++){
        for(i=0;i<t->h;i++){
            for(j=0;j<t->w;j++){

                t->data[l][i][j] = random(mean,var); 
            }
        }
    }

    calcola_stat(t);

}




/**
 * 
 *Matematiche aggiungi e moltiplica scalare
 *
 **/

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
    unsigned int i,j,l;

    ip_mat* output;

    output=ip_mat_create(a->h,a->w,a->k,c);
    
    for(l=0;l<output->k;l++){
        
        for(i=0;i<output->h;i++){

            for(j=0;j<output->w;j++)
            {
                output->data[l][i][j]= output->data[l][i][j] * a->data[l][i][j]; 
            }
        }

    }

    /**
     * aggiorno le statistiche
     * 
     **/
    calcola_stat(output);

    return output;

}


/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output. */
ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
        unsigned int i,j,l;

    ip_mat* output;

    output=ip_mat_create(a->h,a->w,a->k,c);
    
    for(l=0;l<output->k;l++){
        
        for(i=0;i<output->h;i++){

            for(j=0;j<output->w;j++)
            {
                output->data[l][i][j]= output->data[l][i][j] + a->data[l][i][j]; 
            }
        }

    }

    /**
     * aggiorno le statistiche
     * 
     **/
    calcola_stat(output);

    return output;
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

/**
 * Cambiata la funzione get_val: la nostra matrice a indici frastagliati ha l'ordine:
 *      Canale  ->    Colonne   ->  Righe
 * si è modificato solamente l'ordine delle dimensioni
 * 
 * */


float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(k<a->k && i<a->h && j<a->w){  
        return a->data[k][i]][j];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}




/**
float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}
**/

/**
 * Cambiata la funzione set_val: la nostra matrice a indici frastagliati ha l'ordine:
 *      Canale  ->    Colonne   ->  Righe
 * si è modificato solamente l'ordine delle dimensioni
 * 
 * */

/**
void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}
*/
void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(k<a->k && i<a->h && j<a->w){
        a->data[k][i][j]=v;
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
