/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

/*PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT */


/*funzione aux min che trova il minimo di un determinato canale k, riceve in input una ip_mat a e un canale k e restituisce un float che è il minimo  */
float min(ip_mat *a, int k){
    if(a){
        int i,j;
        float minimo;
        minimo = a->data[k][0][0]; /*metto come minimo il primo elemento della matrice*/

        for (i=0; i<a->h; i++){
            for(j=0; j<a->w; j++){ 
                if(minimo > a->data[k][i][j])
                    minimo = a->data[k][i][j];
            }
        }
        return minimo;
    }
    else 
        return 0;
}

/*funzione aux che mi calcola il valore massimo di un determinato canale k*/
float max(ip_mat *a, int k){
    /*verifico se la ip_mat a è diversa da NULL*/
    if(a){
        int i,j;
        float massimo;
        /*all'inizio metto come massimo il primo elemento della matrice e dopo vado a verificare se c'è un altro elemento maggiore*/
        massimo = a->data[k][0][0]; 

        for (i=0; i<a->h; i++){
            for(j=0; j<a->w; j++){
                if(massimo < a->data[k][i][j])
                    massimo = a->data[k][i][j];
            }
        }
        return massimo;
    }
    else 
        return 0;
}

/*mean è la funzione aux che mi permette di calcolare la media degli elementi che si trovano in un determinato canale k*/
float mean(ip_mat *a, int k){
    if(a){ /*verifico se la ip_mat a != NULL*/
        int i,j,nr_el;
        float somma,media;
        somma = 0.0;
        nr_el=0;

        for (i=0; i<a->h; i++){
            for(j=0; j<a->w; j++){
                somma += a->data[k][i][j];
                nr_el++;
            }
        }

        media = somma / nr_el;
        return media;
    }
    else 
        return 0;
}

/* Calcola il valore minimo, il massimo e la media per ogni canale
 * e li salva dentro la struttura ip_mat stats
 * */
void compute_stats(ip_mat * t){
    int i;
    for(i=0; i<t->k; i++){ 
    /*per ogni canale chiamo le funzioni min, max e mean per trovare rispettivamente il minimo, massimo e la media della matrice di quel canale*/

        t->stat->max=max(t,i);
        t->stat->min=min(t,i);
        t->stat->mean=mean(t,i);
    }
}

/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output. */
ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
    if(a->h == b->h && a->w == b->w && a->k == b->k){ /*faccio la verifica per vedere se sono uguali le dimensioni delle ip_mat *a e ip_mat *b */
        int i,j,z;
        ip_mat *out;
        /* siccome ho in input 2 matrici 3D che hanno la stessa dimensione , sommandole ottentiamo una nuova matrice 3D con le stesse dimensioni*/
        out = ip_mat_create(a->h,a->w,a->k,0.0); /*creo la nuova matrice 3D*/
        for (i= 0; i<a->k; i++){
            for (j = 0; j<a->h; j++){
                for (z=0; i<a->w; z++){

                    out-> data[i][j][z] = a->data[i][j][z] + b->data[i][j][z];
                }
            }
        }
        
        compute_stats(out);/*modifico le statistiche per ogni canale della nuova matrice 3D*/

        return out;
    }
    else /*in ramo  else mi trovo sse le due matrici 3D in input hanno dimensione diverse , in questo caso non posso sommare le due matrici a 3 dimensioni*/
        return NULL;
}

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output.
 * */
ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    if(a->h == b->h && a->w == b->w && a->k == b->k){ /* faccio la verifica per vedere se sono uguali le dimensioni delle due matrici a 3 dimensioni */
        int i,j,z;
        ip_mat *out;
        /* siccome ho in input 2 matrici 3D che hanno la stessa dimensione , faccendo la sottrazione ottentiamo una nuova matrice 3D con le stesse dimensioni*/
        out = ip_mat_create(a->h,a->w,a->k,0.0); /*creo la nuova matrice 3D*/
        for (i= 0; i<a->k; i++){
            for (j = 0; j<a->h; j++){
                for (z=0; i<a->w; z++){

                    out-> data[i][j][z] = a->data[i][j][z] - b->data[i][j][z]; /* effetuo la sottrazione */
                }
            }
        }

        compute_stats(out); /*modifico le statistiche per ogni canale della nuova matrice 3D*/

        return out;
    }
    else /*in ramo  else mi trovo sse le due matrici 3D in input hanno dimensione diverse*/
        return NULL;
}

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
    int i,j,z;
    ip_mat *out;
    out = ip_mat_create(a->h,a->w,a->k,0.0); /*creo la nuova matrice 3D, che inizialmente ha tutti i valori 0.0*/

    for (i= 0; i<a->k; i++){
            for (j = 0; j<a->h; j++){
                for (z=0; i<a->w; z++){
                    out-> data[i][j][z] = a->data[i][j][z] * c;
                }
            }
    }

    compute_stats(out); /*modifico le statistiche per ogni canale*/
    return out;
}

/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output. */
ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
    int i,j,z;
    ip_mat *out;
    out = ip_mat_create(a->h,a->w,a->k,0.0); /*creo la nuova matrice 3D*/
    for (i= 0; i<a->k; i++){
            for (j = 0; j<a->h; j++){
                for (z=0; i<a->w; z++){

                    out-> data[i][j][z] = a->data[i][j][z] + c;
                }
            }
    }

    compute_stats(out); /*modifico le statistiche per ogni canale*/

    return out;
}

/* Calcola la media di due ip_mat a e b e la restituisce in output.*/
ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    if(a->h == b->h && a->w == b->w && a->k == b->k){ /*faccio la verifica per vedere se sono uguali le dimensioni delle ip_mat *a e ip_mat *b */
        int i,j,z;
        ip_mat *out;
        /* siccome ho in input 2 matrici 3D che hanno la stessa dimensione , faccendo la media delle 2 , otteniamo una nuova matrice a 3 dimensione che ha le stesse dimensioni*/
        out = ip_mat_create(a->h,a->w,a->k,0.0); /*creo la nuova matrice 3D*/
        for (i= 0; i<a->k; i++){
            for (j = 0; j<a->h; j++){
                for (z=0; i<a->w; z++){

                    out->data[i][j][z] = (a->data[i][j][z] + b->data[i][j][z])/2;
                }
            }
        }

        compute_stats(out); /*modifico le statistiche per ogni canale della nuova matrice a tre dimensioni*/

        return out;
    }
    else /*sse le ip_mat a e b non hanno le stesse dimensioni*/
        return NULL;

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
