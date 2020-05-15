/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

/* HELPERS */
/* le funzioni definite in questa sezione non sono pensate per essere esportate nella libreria, ma sono definite solamente per uso interno */

/* tipo canale: matrice di float */
typedef float **channel_t;
/*
 * posiziona il contenuto del canale source nel canale dest, posizionando la cella 0,0 del canale source alla posizione row,col del canale dest
 * copia il contenuto di source fino a raggiungere il limite del canale dest, e soltanto se row,col sono indici validi per il canale dest
 */
void channel_puts(channel_t dest, unsigned int dest_h, unsigned int dest_w, const channel_t source, unsigned int source_h, unsigned int source_w, unsigned int row, unsigned int col)
{
    if (row < dest_h && col < dest_w)
    {
        unsigned int r, c;
        for (r = 0; r < source_h && r + row < dest_h; r++)
            for (c = 0; c < source_w && c + col < dest_w; c++){
                dest[row + r][col + c] = source[r][c];
            }
    }
}

/*
 * posiziona il contenuto della matrice source nella matrice dest, partendo dalla posizione specificata
 * per tutti i canali di source (massimo il numero di canali di dest) posiziona la cella 0,0 del canale di source alla posizione row,col del canale di dest
 * copia il contenuto di source fino a raggiungere il limite del canale di dest, e soltanto se row,col sono indici validi per la matrice dest
 */
void ip_mat_puts(ip_mat *dest, const ip_mat *source, unsigned int row, unsigned int col)
{
    if (row < dest->h && col < dest->w)
    {
        /* ripeti per tutti i canali di source (massimo dest->k) */
        unsigned int k;
        for (k = 0; k < dest->k && k < source->k; k++)
            channel_puts(dest->data[k], dest->h, dest->w, source->data[k], source->h, source->w, row, col);
    }
}

float easy_random(float mean, float var)
{
    return (mean + var * get_normal_random());
}

/*funzione aux min che trova il minimo di un determinato canale k, riceve in input una ip_mat a e un canale k e restituisce un float che è il minimo  */
float min(ip_mat *a, int k)
{
    if (a)
    {
        unsigned int i, j;
        float minimo;
        minimo = a->data[k][0][0]; /*metto come minimo il primo elemento della matrice*/

        for (i = 0; i < a->h; i++)
        {
            for (j = 0; j < a->w; j++)
            {
                if (minimo > a->data[k][i][j])
                    minimo = a->data[k][i][j];
            }
        }
        return minimo;
    }
    else
        return 0;
}

/*funzione aux che mi calcola il valore massimo di un determinato canale k*/
float max(ip_mat *a, int k)
{
    /*verifico se la ip_mat a è diversa da NULL*/
    if (a)
    {
        unsigned int i, j;
        float massimo;
        /*all'inizio metto come massimo il primo elemento della matrice e dopo vado a verificare se c'è un altro elemento maggiore*/
        massimo = a->data[k][0][0];

        for (i = 0; i < a->h; i++)
        {
            for (j = 0; j < a->w; j++)
            {
                if (massimo < a->data[k][i][j])
                    massimo = a->data[k][i][j];
            }
        }
        return massimo;
    }
    else
        return 0;
}

/*mean è la funzione aux che mi permette di calcolare la media degli elementi che si trovano in un determinato canale k*/
float mean(ip_mat *a, int k)
{
    if (a)
    { /*verifico se la ip_mat a != NULL*/
        unsigned int i, j, nr_el;
        float somma, media;
        somma = 0.0;
        nr_el = 0;

        for (i = 0; i < a->h; i++)
        {
            for (j = 0; j < a->w; j++)
            {
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

/* END HELPERS */

/**** PARTE 1: TIPO DI DATI ip_mat E MEMORIA ****/

/* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */
ip_mat *ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v)
{
    unsigned int i, j, c;

    ip_mat *nuova;

    nuova = (ip_mat *)malloc(sizeof(ip_mat));
    /**
     * asseggno alla struttura nuova i valori dei campi
     * dell'altezza e larghezzza
     *
     **/
    nuova->h = h;
    nuova->k = k;
    nuova->w = w;

    /**creazione delle k matrici + inizializzazione a v**/

    if (nuova)
    {
        /**
         * Creo la matrice a 3 dimensioni effettiva e i 3 vettori con
         * le statistiche per canale
         **/
        nuova->data = (float ***)malloc(sizeof(float **) * k);

        nuova->stat = (stats *)malloc(sizeof(stats) * k);

        /**
             * canale
             **/
        for (c = 0; c < k; c++)
        {

            nuova->data[c] = (float **)malloc(sizeof(float *) * h);
        }

        for (c = 0; c < k; c++)
        {

            for (i = 0; i < h; i++)
            {

                nuova->data[c][i] = (float *)malloc(sizeof(float) * w);
            }
        }

        /*inzializzazione dei valori della matrice*/

        for (c = 0; c < k; c++)
        {
            for (i = 0; i < h; i++)
            {
                for (j = 0; j < w; j++)
                {
                    nuova->data[c][i][j] = v;
                }
            }
        }

        /**
             * riempio i valri di stats per ogni canale;
             * 
             **/

        for (c = 0; c < k; c++)
        {
            nuova->stat[c].max = v;
            nuova->stat[c].min = v;
            nuova->stat[c].mean = mean(nuova, (int)c);
        }
    }

    return nuova;
}

/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat *a)
{
    if (a != NULL)
    {
        unsigned int ch, row;
        /* libera stats */
        free(a->stat);
        /* libera canali della matrice data */
        for (ch = 0; ch < a->k; ch++)
        {
            /* libera ogni riga del canale: array di colonne */
            for (row = 0; row < a->h; row++)
                free(a->data[ch][row]);
            /* libera il canale: array di righe */
            free(a->data[ch]);
        }
        /* libera data */
        free(a->data);
        /* libera la struttura */
        free(a);
    }
}

/**
 * Cambiata la funzione get_val: la nostra matrice a indici frastagliati ha l'ordine:
 *      Canale  ->    Colonne   ->  Righe
 * si è modificato solamente l'ordine delle dimensioni
 * 
 * */

float get_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k)
{
    if (k < a->k && i < a->h && j < a->w)
    {
        return a->data[k][i][j];
    }
    else
    {
        printf("Errore get_val!!!");
        exit(1);
    }
}

/**
 * Cambiata la funzione set_val: la nostra matrice a indici frastagliati ha l'ordine:
 *      Canale  ->    Colonne   ->  Righe
 * si è modificato solamente l'ordine delle dimensioni
 * 
 * */

void set_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k, float v)
{
    if (k < a->k && i < a->h && j < a->w)
    {
        a->data[k][i][j] = v;
    }
    else
    {
        printf("Errore set_val!!!");
        exit(1);
    }
}

/* Calcola il valore minimo, il massimo e la media per ogni canale
 * e li salva dentro la struttura ip_mat stats
 * */
void compute_stats(ip_mat *t)
{
    unsigned int i;
    for (i = 0; i < t->k; i++)
    {
        /*per ogni canale chiamo le funzioni min, max e mean per trovare rispettivamente il minimo, massimo e la media della matrice di quel canale*/

        t->stat->max = max(t, i);
        t->stat->min = min(t, i);
        t->stat->mean = mean(t, i);
    }
}

/**
 * TODO: random(mean, var); 
 * Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e varianza var 
 * 
 * */

void ip_mat_init_random(ip_mat *t, float mean, float var)
{
    unsigned int i, j, l;

    /**
     * 
     * riempio la matrice del canale l-esimo con valori casuali
     *
     * */

    for (l = 0; l < t->k; l++)
    {
        for (i = 0; i < t->h; i++)
        {
            for (j = 0; j < t->w; j++)
            {

                t->data[l][i][j] = easy_random(mean, var);
            }
        }
    }

    compute_stats(t);
}

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat *ip_mat_copy(ip_mat *in)
{
    return ip_mat_subset(in, 0, in->h, 0, in->w);
}

/* Restituisce una sotto-matrice, ovvero la porzione individuata da:
 * t->data[row_start...row_end][col_start...col_end][0...k]
 * La terza dimensione la riportiamo per intero, stiamo in sostanza prendendo un sottoinsieme
 * delle righe e delle colonne.
 * */
ip_mat *ip_mat_subset(ip_mat *t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end)
{
    ip_mat *subset_mat = NULL;
    if (row_start <= row_end && col_start <= col_end && row_end <= t->h && col_end <= t->w)
    {
        unsigned int ch, row, col;
        subset_mat = ip_mat_create(row_end - row_start, col_end - col_start, t->k, 0.0);
        for (ch = 0; ch < t->k; ch++)
        {
            /* copia gli stats per il canale */
            subset_mat->stat[ch] = t->stat[ch];
            /* copia i dati per il canale */
            for (row = row_start; row < row_end; row++)
                for (col = col_start; col < col_end; col++)
                    subset_mat->data[ch][row][col] = t->data[ch][row][col];
        }
    }
    return subset_mat;
}

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
 * */
ip_mat *ip_mat_concat(ip_mat *a, ip_mat *b, int dimensione)
{
    ip_mat *concat_mat = NULL;
    switch (dimensione)
    {
    case 0:
        /* altezza */
        if (a->w == b->w && a->k == b->k)
        {
            concat_mat = ip_mat_create(a->h + b->h, a->w, a->k, 0.0);
            if (concat_mat)
            {
                ip_mat_puts(concat_mat, a, 0, 0);
                ip_mat_puts(concat_mat, b, a->h, 0);
            }
        }
        break;
    case 1:
        /* larghezza */
        if (a->h == b->h && a->k == b->k)
        {
            concat_mat = ip_mat_create(a->h, a->w + b->w, a->k, 0.0);
            if (concat_mat)
            {
                ip_mat_puts(concat_mat, a, 0, 0);
                ip_mat_puts(concat_mat, b, 0, a->w);
            }
        }
        break;
    case 2:
        /* canale */
        if (a->w == b->w && a->h == b->h)
        {
            concat_mat = ip_mat_create(a->h, a->w, a->k + b->k, 0.0);
            if (concat_mat)
            {
                /* posiziona i canali di a */
                ip_mat_puts(concat_mat, a, 0, 0);
                /* sposta temporaneamente in avanti il puntatore dei canali al canale k, per permettere di caricare i canali di b */
                concat_mat->data += a->k;
                /* posiziona i canali di b */
                ip_mat_puts(concat_mat, b, 0, 0);
                /* riporta il puntatore dei canali al canale 0 */
                concat_mat->data -= a->k;
            }
        }
        break;
    default:
        break;
    }
    return concat_mat;
}

/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/

/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output. */
ip_mat *ip_mat_sum(ip_mat *a, ip_mat *b)
{
    if (a->h == b->h && a->w == b->w && a->k == b->k)
    { /*faccio la verifica per vedere se sono uguali le dimensioni delle ip_mat *a e ip_mat *b */
        unsigned int i, j, z;
        ip_mat *out;
        /* siccome ho in input 2 matrici 3D che hanno la stessa dimensione , sommandole ottentiamo una nuova matrice 3D con le stesse dimensioni*/
        out = ip_mat_create(a->h, a->w, a->k, 0.0); /*creo la nuova matrice 3D*/
        for (i = 0; i < a->k; i++)
        {
            for (j = 0; j < a->h; j++)
            {
                for (z = 0; i < a->w; z++)
                {

                    out->data[i][j][z] = a->data[i][j][z] + b->data[i][j][z];
                }
            }
        }

        compute_stats(out); /*modifico le statistiche per ogni canale della nuova matrice 3D*/

        return out;
    }
    else /*in ramo  else mi trovo sse le due matrici 3D in input hanno dimensione diverse , in questo caso non posso sommare le due matrici a 3 dimensioni*/
        return NULL;
}

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output.
 * */
ip_mat *ip_mat_sub(ip_mat *a, ip_mat *b)
{
    if (a->h == b->h && a->w == b->w && a->k == b->k)
    { /* faccio la verifica per vedere se sono uguali le dimensioni delle due matrici a 3 dimensioni */
        unsigned int i, j, z;
        ip_mat *out;
        /* siccome ho in input 2 matrici 3D che hanno la stessa dimensione , faccendo la sottrazione ottentiamo una nuova matrice 3D con le stesse dimensioni*/
        out = ip_mat_create(a->h, a->w, a->k, 0.0); /*creo la nuova matrice 3D*/
        for (i = 0; i < a->k; i++)
        {
            for (j = 0; j < a->h; j++)
            {
                for (z = 0; i < a->w; z++)
                {

                    out->data[i][j][z] = a->data[i][j][z] - b->data[i][j][z]; /* effetuo la sottrazione */
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
ip_mat *ip_mat_mul_scalar(ip_mat *a, float c)
{
    unsigned int i, j, z;
    ip_mat *out;
    out = ip_mat_create(a->h, a->w, a->k, 0.0); /*creo la nuova matrice 3D, che inizialmente ha tutti i valori 0.0*/

    for (i = 0; i < a->k; i++)
    {
        for (j = 0; j < a->h; j++)
        {
            for (z = 0; i < a->w; z++)
            {
                out->data[i][j][z] = a->data[i][j][z] * c;
            }
        }
    }

    compute_stats(out); /*modifico le statistiche per ogni canale*/
    return out;
}

/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output. */
ip_mat *ip_mat_add_scalar(ip_mat *a, float c)
{
    unsigned int ch, row, col;
    ip_mat *out;
    out = ip_mat_create(a->h, a->w, a->k, 0.0); /*creo la nuova matrice 3D*/
    for (ch = 0; ch < a->k; ch++)
    {
        for (row = 0; row < a->h; row++)
        {
            for (col = 0; col < a->w; col++)
            {

                out->data[ch][row][col] = a->data[ch][row][col] + c;
            }
        }
    }

    compute_stats(out); /*modifico le statistiche per ogni canale*/

    return out;
}

/* Calcola la media di due ip_mat a e b e la restituisce in output.*/
ip_mat *ip_mat_mean(ip_mat *a, ip_mat *b)
{
    if (a->h == b->h && a->w == b->w && a->k == b->k)
    { /*faccio la verifica per vedere se sono uguali le dimensioni delle ip_mat *a e ip_mat *b */
        unsigned int i, j, z;
        ip_mat *out;
        /* siccome ho in input 2 matrici 3D che hanno la stessa dimensione , faccendo la media delle 2 , otteniamo una nuova matrice a 3 dimensione che ha le stesse dimensioni*/
        out = ip_mat_create(a->h, a->w, a->k, 0.0); /*creo la nuova matrice 3D*/
        for (i = 0; i < a->k; i++)
        {
            for (j = 0; j < a->h; j++)
            {
                for (z = 0; i < a->w; z++)
                {

                    out->data[i][j][z] = (a->data[i][j][z] + b->data[i][j][z]) / 2;
                }
            }
        }

        compute_stats(out); /*modifico le statistiche per ogni canale della nuova matrice a tre dimensioni*/

        return out;
    }
    else /*sse le ip_mat a e b non hanno le stesse dimensioni*/
        return NULL;
}

/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/
/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 * */
ip_mat *ip_mat_to_gray_scale(ip_mat *in)
{
    return NULL;
}

/* Effettua la fusione (combinazione convessa) di due immagini */
ip_mat *ip_mat_blend(ip_mat *a, ip_mat *b, float alpha)
{
    return NULL;
}

/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore*/
ip_mat *ip_mat_brighten(ip_mat *a, float bright)
{
    return ip_mat_add_scalar(a, bright);
}

/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 * out = a + gauss_noise*amount
 * */
ip_mat *ip_mat_corrupt(ip_mat *a, float amount)
{
    return NULL;
}

/**** PARTE 3: CONVOLUZIONE E FILTRI *****/

/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a".
 * */
ip_mat *ip_mat_convolve(ip_mat *a, ip_mat *f)
{
    return NULL;
}

/* Aggiunge un padding all'immagine. Il padding verticale è pad_h mentre quello
 * orizzontale è pad_w.
 * L'output sarà un'immagine di dimensioni:
 *      out.h = a.h + 2*pad_h;
 *      out.w = a.w + 2*pad_w;
 *      out.k = a.k
 * con valori nulli sui bordi corrispondenti al padding e l'immagine "a" riportata
 * nel centro
 * */
ip_mat *ip_mat_padding(ip_mat *a, int pad_h, int pad_w)
{
    return NULL;
}

/* Crea un filtro di sharpening */
ip_mat *create_sharpen_filter()
{
    return NULL;
}

/* Crea un filtro per rilevare i bordi */
ip_mat *create_edge_filter()
{
    return NULL;
}

/* Crea un filtro per aggiungere profondità */
ip_mat *create_emboss_filter()
{
    return NULL;
}

/* Crea un filtro medio per la rimozione del rumore */
ip_mat *create_average_filter(int w, int h, int k)
{
    return NULL;
}

/* Crea un filtro gaussiano per la rimozione del rumore */
ip_mat *create_gaussian_filter(int w, int h, int k, float sigma)
{
    return NULL;
}

/* Effettua una riscalatura dei dati tale che i valori siano in [0,new_max].
 * Utilizzate il metodo compute_stat per ricavarvi il min, max per ogni canale.
 *
 * I valori sono scalati tramite la formula valore-min/(max - min)
 *
 * Si considera ogni indice della terza dimensione indipendente, quindi l'operazione
 * di scalatura va ripetuta per ogni "fetta" della matrice 3D.
 * Successivamente moltiplichiamo per new_max gli elementi della matrice in modo da ottenere un range
 * di valori in [0,new_max].
 * */
void rescale(ip_mat *t, float new_max)
{
}

/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.*/
void clamp(ip_mat *t, float low, float high)
{
}

/**** METODI GIA' IMPLEMENTATI ****/

/* Genera dei numeri casuali con distribuzione Normale (versione base)
 * https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
 * */
float get_normal_random()
{
    float y1 = ((float)(rand()) + 1.) / ((float)(RAND_MAX) + 1.);
    float y2 = ((float)(rand()) + 1.) / ((float)(RAND_MAX) + 1.);
    return cos(2 * PI * y2) * sqrt(-2. * log(y1));
}

/* Converte una Bitmap in una ip_mat*/
ip_mat *bitmap_to_ip_mat(Bitmap *img)
{
    unsigned int i = 0, j = 0;

    unsigned char R, G, B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat *out = ip_mat_create(h, w, 3, 0);

    for (i = 0; i < h; i++)
    {
        for (j = 0; j < w; j++)
        {
            bm_get_pixel(img, j, i, &R, &G, &B);
            set_val(out, i, j, 0, (float)R);
            set_val(out, i, j, 1, (float)G);
            set_val(out, i, j, 2, (float)B);
        }
    }

    return out;
}

/* Converte una ip_mat in una bitmap*/
Bitmap *ip_mat_to_bitmap(ip_mat *t)
{

    Bitmap *b = bm_create(t->w, t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)
    {
        for (j = 0; j < t->w; j++)
        {
            bm_set_pixel(b, j, i, (unsigned char)get_val(t, i, j, 0),
                         (unsigned char)get_val(t, i, j, 1),
                         (unsigned char)get_val(t, i, j, 2));
        }
    }
    return b;
}

/* Visualizza i dati stampando in ordine le matrici rispetto
 * la terza dimensione.
 * Prima stamperemo t->data[...][...][0] poi t->data[...][...][1] ...
 * */
void ip_mat_show(ip_mat *t)
{
    unsigned int i, l, j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n", t->h, t->w, t->k);
    for (l = 0; l < t->k; l++)
    {
        printf("Slice %d\n", l);
        for (i = 0; i < t->h; i++)
        {
            for (j = 0; j < t->w; j++)
            {
                printf("%f ", get_val(t, i, j, l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

/* Visualizza a video le statistiche per ogni canale.
 * */
void ip_mat_show_stats(ip_mat *t)
{
    unsigned int k;

    for (k = 0; k < t->k; k++)
    {
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}