/**
*  CIMP - C Image Manipulation Program
*  (https://github.com/EPMatt/cimp)
*
*  (C) 2020
*
*  Ina Popescu   (https://github.com/ina-pps)
*  Matteo Agnoletto   (https://github.com/EPMatt)
*  Lorenzo Donatelli  (https://github.com/whitedemond)
*  
*  For licensing conditions related to this project, see LICENSE
*
*/

/*
 Laboratorio di programmazione A.A. 2019/2020

 Sebastiano Vascon
*/

#ifndef IP_LIB_H
#define IP_LIB_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "bmp.h"

#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */
#define PI 3.141592654

typedef struct
{
    float min;
    float max;
    float mean;
} stats;

typedef struct
{
    unsigned int w; /* <- larghezza */
    unsigned int h; /* <- altezza */
    unsigned int k; /* <- canali */
    stats *stat;    /* <- statistiche per canale */
    float ***data;  /* <- matrice 3D di valori float */
} ip_mat;

/* HELPERS */
/* le funzioni, i tipi e le costanti sono pensate per utilizzo interno alla libreria */

/* valori massimi e minimi, in virgola mobile, che un pixel può assumere */
#define MAX_PIXEL_FLOAT 255.0
#define MIN_PIXEL_FLOAT 0.0

/* costanti utilizzati nel richiamare la funzione ip_mat_puts */
#define NO_COMPUTE_STATS 0
#define COMPUTE_STATS 1

/* tipo canale: matrice di float */
typedef struct
{
    float **data;
    unsigned int h;
    unsigned int w;
} channel_t;

/* data una ip_mat, ottieni una struttura canale, con riferimento al ch-esimo canale */
channel_t get_channel(ip_mat const *a, unsigned int ch);

/*
 * posiziona il contenuto del canale source nel canale dest, posizionando la cella 0,0 del canale source alla posizione row,col del canale dest
 * copia il contenuto di source fino a raggiungere il limite del canale dest, e soltanto se row,col sono indici validi per il canale dest
 * */
void channel_puts(channel_t const dest, channel_t const source, unsigned int row, unsigned int col);

/*
 * posiziona il contenuto della matrice source nella matrice dest, partendo dalla posizione specificata
 * per tutti i canali di source (massimo il numero di canali di dest) posiziona la cella 0,0 del canale di source alla posizione row,col del canale di dest
 * copia il contenuto di source fino a raggiungere il limite del canale di dest, e soltanto se row,col sono indici validi per la matrice dest
 * */
void ip_mat_puts(ip_mat *dest, ip_mat const *source, unsigned int row, unsigned int col, int do_compute_stats);

/**
 * restituisce un messaggio di errore nel caso in cui la matrice passata sia NULL
 **/
void not_null_ip_mat(ip_mat const *a);

/*funzione aux min che trova il minimo di un determinato canale k, riceve in input una ip_mat a e un canale k e restituisce un float che è il minimo  */
float min(ip_mat *a, unsigned int k);

/*funzione aux che mi calcola il valore massimo di un determinato canale k*/
float max(ip_mat *a, unsigned int k);

/*mean è la funzione aux che mi permette di calcolare la media degli elementi che si trovano in un determinato canale k*/
float mean(ip_mat *a, unsigned int k);

/* restringi il valore fornito all'interno del range low...high, estremi inclusi */
float restrict_val(float val, float low, float high);

/**
 * Calcola la media per il pixel (i,j) sui k canali della matrice passata come parametro.
 * 
 * */
float mean_pixel_channel(ip_mat *a, unsigned int i, unsigned int j);

/**
 * Controlla se due ip_mat hanno le stesse dimensioni oppure no
 * esce dal programma nel momento in cui una delle dimensioni é diversa
 */
void equal_dimension(ip_mat *a, ip_mat *b);

/* calcola la somma di prodotti tra il kernel fornito e il canale fornito, partendo dalla posizione (start_h,start_w) del canale */
float convolve_channel(channel_t ch, channel_t filter, unsigned int start_h, unsigned int start_w);

/* END HELPERS */

/**** PARTE 1: TIPO DI DATI ip_mat E MEMORIA ****/

/* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */
ip_mat *ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v);

/* Libera la memoria (data, stat e la struttura)
 *
 * se la variabile "a" è NULL non fa nulla.
 *
 * */
void ip_mat_free(ip_mat *a);

/* Restituisce il valore in posizione i,j,k */
float get_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k);

/* Setta il valore in posizione i,j,k a v*/
void set_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k, float v);

/* Calcola il valore minimo, il massimo e la media per ogni canale
 * e li salva dentro la struttura ip_mat stats
 * */
void compute_stats(ip_mat *t);

/* Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e deviazione std */
void ip_mat_init_random(ip_mat *t, float mean, float std);

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat *ip_mat_copy(ip_mat *in);

/* Restituisce una sotto-matrice, ovvero la porzione individuata da:
 * t->data[row_start...row_end][col_start...col_end][0...k]
 * La terza dimensione la riportiamo per intero, stiamo in sostanza prendendo un sottoinsieme
 * delle righe e delle colonne.
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 * */
ip_mat *ip_mat_subset(ip_mat *t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end);

/* Concatena due ip_mat su una certa dimensione.
 * Ad esempio:
 * ip_mat_concat(ip_mat * a, ip_mat * b, 0);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h + b.h
 *      out.w = a.w = b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 1);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w + b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 2);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w = b.w
 *      out.k = a.k + b.k
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 * */
ip_mat *ip_mat_concat(ip_mat *a, ip_mat *b, int dimensione);

/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/
/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche).
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_sum(ip_mat *a, ip_mat *b);

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_sub(ip_mat *a, ip_mat *b);

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_mul_scalar(ip_mat *a, float c);

/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output.
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_add_scalar(ip_mat *a, float c);

/* Calcola la media di due ip_mat a e b. La media si calcola per coppie delle due matrici aventi gli stessi indici
 * C[i][j][k]= (A[i][j][k]+B[i]j[k])/2
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_mean(ip_mat *a, ip_mat *b);

/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/
/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 * */
ip_mat *ip_mat_to_gray_scale(ip_mat *in);

/* Effettua la fusione (combinazione convessa) di due immagini.
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 *
 * Le variabili "a" e "b" devono avere le stesse dimensioni
 */
ip_mat *ip_mat_blend(ip_mat *a, ip_mat *b, float alpha);

/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_brighten(ip_mat *a, float bright);

/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 *
 * out = a + gauss_noise*amount
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 * */
ip_mat *ip_mat_corrupt(ip_mat *a, float amount);

/**** PARTE 3: CONVOLUZIONE E FILTRI *****/

/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a".
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_convolve(ip_mat *a, ip_mat *f);

/* Aggiunge un padding all'immagine. Il padding verticale è pad_h mentre quello
 * orizzontale è pad_w.
 * L'output sarà un'immagine di dimensioni:
 *      out.h = a.h + 2*pad_h;
 *      out.w = a.w + 2*pad_w;
 *      out.k = a.k
 * con valori nulli sui bordi corrispondenti al padding e l'immagine "a" riportata
 * nel centro.
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_padding(ip_mat *a, unsigned int pad_h, unsigned int pad_w);

/* Crea un filtro di sharpening */
ip_mat *create_sharpen_filter();

/* Crea un filtro per rilevare i bordi */
ip_mat *create_edge_filter();

/* Crea un filtro per aggiungere profondità */
ip_mat *create_emboss_filter();

/* Crea un filtro medio per la rimozione del rumore */
ip_mat *create_average_filter(unsigned int h, unsigned int w, unsigned int k);

/* Crea un filtro gaussiano per la rimozione del rumore */
ip_mat *create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma);

/* Effettua una riscalatura dei dati tale che i valori siano in [0,new_max].
 * Utilizzate il metodo compute_stat per ricavarvi il min, max per ogni canale.
 *
 * I valori sono scalati tramite la formula (valore-min)/(max - min)
 *
 * Si considera ogni indice della terza dimensione indipendente, quindi l'operazione
 * di scalatura va ripetuta per ogni "fetta" della matrice 3D.
 * Successivamente moltiplichiamo per new_max gli elementi della matrice in modo da ottenere un range
 * di valori in [0,new_max].
 *
 * Il risultato dell'operazione si salva in t
 * */
void rescale(ip_mat *t, float new_max);

/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.
 *
 * Il risultato dell'operazione si salva in t
 * */
void clamp(ip_mat *t, float low, float high);

/**** METODI GIA' IMPLEMENTATI ****/
/* Genera dei numeri casuali con distribuzione Normale
 * https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
 * caratterizzata da media = mu e deviazione  = std
 * */
float get_normal_random(float media, float std);

/* Converte una Bitmap in una ip_mat*/
ip_mat *bitmap_to_ip_mat(Bitmap *img);

/* Converte una ip_mat in una bitmap*/
Bitmap *ip_mat_to_bitmap(ip_mat *t);

/* Visualizza i dati stampando in ordine le matrici rispetto
 * la terza dimensione.
 * Prima stamperemo t->data[...][...][0] poi t->data[...][...][1] ...
 * */
void ip_mat_show(ip_mat *t);

/* Visualizza a video le statistiche per ogni canale.
 * */
void ip_mat_show_stats(ip_mat *t);

#endif /*IP_LIB_H*/
