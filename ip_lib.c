/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

/* HELPERS */

/* data una ip_mat, ottieni una struttura canale, con riferimento al ch-esimo canale */
channel_t get_channel(ip_mat const *a, unsigned int ch)
{
    not_null_ip_mat(a);
    if (ch < a->k)
    {
        channel_t out;
        out.h = a->h;
        out.w = a->w;
        /* puntatore al canale */
        out.data = a->data[ch];
        return out;
    }
    else
    {
        printf("Invalid channel number\n");
        exit(1);
    }
}

/*
 * posiziona il contenuto del canale source nel canale dest, posizionando la cella 0,0 del canale source alla posizione row,col del canale dest
 * copia il contenuto di source fino a raggiungere il limite del canale dest, e soltanto se row,col sono indici validi per il canale dest
 * */
void channel_puts(channel_t const dest, channel_t const source, unsigned int row, unsigned int col)
{
    if (row < dest.h && col < dest.w)
    {
        unsigned int r, c;
        for (r = 0; r < source.h && r + row < dest.h; r++)
            for (c = 0; c < source.w && c + col < dest.w; c++)
                dest.data[row + r][col + c] = source.data[r][c];
    }
}

/*
 * posiziona il contenuto della matrice source nella matrice dest, partendo dalla posizione specificata
 * per tutti i canali di source (massimo il numero di canali di dest) posiziona la cella 0,0 del canale di source alla posizione row,col del canale di dest
 * copia il contenuto di source fino a raggiungere il limite del canale di dest, e soltanto se row,col sono indici validi per la matrice dest
 * */
void ip_mat_puts(ip_mat *dest, ip_mat const *source, unsigned int row, unsigned int col, int do_compute_stats)
{
    not_null_ip_mat(dest);
    not_null_ip_mat(source);
    if (row < dest->h && col < dest->w)
    {
        /* ripeti per tutti i canali di source (massimo dest->k) */
        unsigned int k;
        for (k = 0; k < dest->k && k < source->k; k++)
            channel_puts(get_channel(dest, k), get_channel(source, k), row, col);
        /* se richiesto ricalcola gli stats */
        if (do_compute_stats)
            compute_stats(dest);
    }
}

/**
 * restituisce un messaggio di errore nel caso in cui la matrice passata sia NULL
 **/
void not_null_ip_mat(ip_mat const *a)
{
    if (a == NULL)
    {
        printf("Parameter passed NULL\n");
        exit(1);
    }
}

/*funzione aux min che trova il minimo di un determinato canale k, riceve in input una ip_mat a e un canale k e restituisce un float che è il minimo  */
float min(ip_mat *a, int k)
{
    unsigned int i, j;
    float minimo = 0.0, current;

    not_null_ip_mat(a);

    minimo = get_val(a, 0, 0, k); /*metto come minimo il primo elemento della matrice*/

    for (i = 0; i < a->h; i++)
    {
        for (j = 0; j < a->w; j++)
        {
            current = get_val(a, i, j, k);
            if (minimo > current)
                minimo = current;
        }
    }
    return minimo;
}

/*funzione aux che mi calcola il valore massimo di un determinato canale k*/
float max(ip_mat *a, int k)
{
    unsigned int i, j;
    float massimo = 0.0, current;

    /*verifico se la ip_mat a è diversa da NULL*/

    not_null_ip_mat(a);

    /*all'inizio metto come massimo il primo elemento della matrice e dopo vado a verificare se c'è un altro elemento maggiore*/

    massimo = get_val(a, 0, 0, k);

    for (i = 0; i < a->h; i++)
    {
        for (j = 0; j < a->w; j++)
        {
            current = get_val(a, i, j, k);
            if (massimo < current)
                massimo = current;
        }
    }
    return massimo;
}

/*mean è la funzione aux che mi permette di calcolare la media degli elementi che si trovano in un determinato canale k*/
float mean(ip_mat *a, int k)
{
    unsigned int row, col;
    float somma = 0.0;
    float media = 0.0;

    /*verifico se la ip_mat a != NULL*/
    not_null_ip_mat(a);

    for (row = 0; row < a->h; row++)
    {
        for (col = 0; col < a->w; col++)
            somma += get_val(a, row, col, k);
    }

    media = somma / (a->h * a->w);
    return media;
}

/* restringi il valore fornito all'interno del range low...high, estremi inclusi */
float restrict_val(float val, float low, float high)
{
    if (val < low)
        return low;
    else if (val > high)
        return high;
    else
        return val;
}

/**
 * Calcola la media per ogni pixel sui k canali della matrice di pixel, prendendo un indice colonna e un indice riga.
 * 
 * */
float mean_pixel_channel(ip_mat *a, unsigned int i, unsigned int j)
{
    unsigned int ch;
    float sup = 0.0;

    not_null_ip_mat(a);

    for (ch = 0; ch < a->k; ch++)
        sup += get_val(a, i, j, ch);

    sup /= a->k;

    return sup;
}

/**
 * Controlla se due ip_mmat hanno le stesse dimensioni oppre no
 * esce dal programma nel momento in cui una delle dimensioni é diversa
 */
void equal_dimension(ip_mat *a, ip_mat *b)
{
    not_null_ip_mat(a);
    not_null_ip_mat(b);

    if ((a->h != b->h) || ((a->w != b->w) && (a->k != b->k)))
    {
        printf("No equal size of images! \nPlease, Provide images with the same dimensions \n");
        exit(1);
    }
}

/**
 * restituisce un messaggio di errore nel caso in cui la matrice passata sia NULL
 **/

/* calcola la somma di prodotti tra il kernel fornito e il canale fornito, partendo dalla posizione (start_h,start_w) del canale */
float convolve_channel(channel_t ch, channel_t filter, unsigned int start_h, unsigned int start_w)
{
    if (start_h + filter.h <= ch.h && start_w + filter.w <= ch.w)
    {
        float result = 0.0;
        unsigned int row, col;
        for (row = 0; row < filter.h; row++)
            for (col = 0; col < filter.w; col++)
                result += ch.data[start_h + row][start_w + col] * filter.data[row][col];
        return result;
    }
    else
    {
        printf("Invalid start_h or start_w \n");
        exit(1);
    }
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
     * assegno alla struttura nuova i valori dei campi
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
        /* canale */

        for (c = 0; c < k; c++)
            nuova->data[c] = (float **)malloc(sizeof(float *) * h);
        for (c = 0; c < k; c++)
        {
            if (nuova->data[c])
            {
                for (i = 0; i < h; i++)
                    nuova->data[c][i] = (float *)malloc(sizeof(float) * w);
            }
            else
            {
                printf("Unable to allocate nuova->data[c]\n c generic position\n");
                exit(1);
            }
        }
        /* inizializzazione dei valori della matrice */
        for (c = 0; c < k; c++)
        {
            for (i = 0; i < h; i++)
            {
                if (nuova->data[c][i])
                {
                    for (j = 0; j < w; j++)
                        set_val(nuova, i, j, c, v);
                }
                else
                {
                    printf("Unable to allocate nuova->data[c][i]\n c,i generic position\n");
                    exit(1);
                }
            }
        }
        /* riempio i valori di stats per ogni canale */
        if (nuova->stat)
        {
            for (c = 0; c < k; c++)
            {
                nuova->stat[c].max = v;
                nuova->stat[c].min = v;
                nuova->stat[c].mean = v;
            }
        }

        return nuova;
    }
    else
    { /** nel caso in cui la creazione della nuova ip_mat fallisca **/

        printf("Unable to allocate matrix structure");
        exit(1);
    }
}

/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat *a)
{
    if (a)
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

/* Restituisce il valore in posizione i,j,k */
float get_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k)
{
    not_null_ip_mat(a);
    if (i < a->h && j < a->w && k < a->k)
        return a->data[k][i][j];
    else
    {
        printf("Errore get_val!!!");
        exit(1);
    }
}

/* Setta il valore in posizione i,j,k a v*/
void set_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k, float v)
{
    not_null_ip_mat(a);
    if (i < a->h && j < a->w && k < a->k)
        a->data[k][i][j] = v;
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
    not_null_ip_mat(t);

    for (i = 0; i < t->k; i++)
    {
        /*per ogni canale chiamo le funzioni min, max e mean per trovare rispettivamente il minimo, massimo e la media della matrice di quel canale*/
        t->stat[i].max = max(t, i);
        t->stat[i].min = min(t, i);
        t->stat[i].mean = mean(t, i);
    }
}

/* Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e deviazione std */
void ip_mat_init_random(ip_mat *t, float mean, float std)
{
    unsigned int row, col, ch;

    not_null_ip_mat(t);
    /**
     * 
     * riempio la matrice del canale l-esimo con valori casuali
     *
     * */
    for (ch = 0; ch < t->k; ch++)
    {
        for (row = 0; row < t->h; row++)
        {
            for (col = 0; col < t->w; col++)
                set_val(t, row, col, ch, get_normal_random(mean, std));
        }
    }
    compute_stats(t);
}

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat *ip_mat_copy(ip_mat *in)
{
    unsigned int ch;
    ip_mat *new_mat;

    not_null_ip_mat(in);

    new_mat = ip_mat_create(in->h, in->w, in->k, 0.0);

    ip_mat_puts(new_mat, in, 0, 0, NO_COMPUTE_STATS);
    /* copia gli stats dei canali senza ricalcolarli */
    for (ch = 0; ch < new_mat->k; ch++)
        new_mat->stat[ch] = in->stat[ch];

    return new_mat;
}

/* Restituisce una sotto-matrice, ovvero la porzione individuata da:
 * t->data[row_start...row_end][col_start...col_end][0...k]
 * La terza dimensione la riportiamo per intero, stiamo in sostanza prendendo un sottoinsieme
 * delle righe e delle colonne.
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 * */
ip_mat *ip_mat_subset(ip_mat *t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end)
{
    ip_mat *subset_mat = NULL;

    not_null_ip_mat(t);

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
                    set_val(subset_mat, row, col, ch, get_val(t, row, col, ch));
        }
    }
    compute_stats(subset_mat);
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
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 * */
ip_mat *ip_mat_concat(ip_mat *a, ip_mat *b, int dimensione)
{
    ip_mat *concat_mat;
    int valid_dim = 1;
    not_null_ip_mat(a);
    not_null_ip_mat(b);
    switch (dimensione)
    {
    case 0:
        /* altezza */
        valid_dim = a->w == b->w && a->k == b->k;
        if (valid_dim)
        {
            concat_mat = ip_mat_create(a->h + b->h, a->w, a->k, 0.0);
            ip_mat_puts(concat_mat, a, 0, 0, NO_COMPUTE_STATS);
            ip_mat_puts(concat_mat, b, a->h, 0, COMPUTE_STATS);
        }
        break;
    case 1:
        /* larghezza */
        valid_dim = a->h == b->h && a->k == b->k;
        if (valid_dim)
        {
            concat_mat = ip_mat_create(a->h, a->w + b->w, a->k, 0.0);
            ip_mat_puts(concat_mat, a, 0, 0, NO_COMPUTE_STATS);
            ip_mat_puts(concat_mat, b, 0, a->w, COMPUTE_STATS);
        }
        break;
    case 2:
        /* canale */
        valid_dim = a->w == b->w && a->h == b->h;
        if (valid_dim)
        {
            concat_mat = ip_mat_create(a->h, a->w, a->k + b->k, 0.0);
            /* posiziona i canali di a */
            ip_mat_puts(concat_mat, a, 0, 0, NO_COMPUTE_STATS);
            /* sposta temporaneamente in avanti il puntatore dei canali al canale k, e riduci il numero di canali, per permettere di caricare i canali di b */
            concat_mat->data += a->k;
            concat_mat->k -= a->k;
            /* posiziona i canali di b */
            ip_mat_puts(concat_mat, b, 0, 0, NO_COMPUTE_STATS);
            /* riporta il puntatore dei canali al canale 0, e il numero di canali a quello completo */
            concat_mat->data -= a->k;
            concat_mat->k += a->k;
            /* ricalcola gli stats su tutti i canali */
            compute_stats(concat_mat);
        }

        break;
    default:
        printf("Invalid concat option\n");
        exit(1);
        break;
    }
    if (valid_dim)
        return concat_mat;
    else
    {
        printf("Invalid dimensions for required operation\n");
        exit(1);
    }
}

/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/
/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche).
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_sum(ip_mat *a, ip_mat *b)
{
    unsigned int ch, row, col;
    ip_mat *out;

    not_null_ip_mat(a);
    not_null_ip_mat(b);

    /*faccio la verifica per vedere se sono uguali le dimensioni delle ip_mat *a e ip_mat *b */
    equal_dimension(a, b);
    /* siccome ho in input 2 matrici 3D che hanno la stessa dimensione , sommandole ottentiamo una nuova matrice 3D con le stesse dimensioni*/
    out = ip_mat_create(a->h, a->w, a->k, 0.0); /*creo la nuova matrice 3D*/
    for (ch = 0; ch < a->k; ch++)
    {
        for (row = 0; row < a->h; row++)
        {
            for (col = 0; col < a->w; col++)
                set_val(out, row, col, ch, (get_val(a, row, col, ch) + get_val(b, row, col, ch)));
        }
    }
    compute_stats(out); /*modifico le statistiche per ogni canale della nuova matrice 3D*/
    return out;
}

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_sub(ip_mat *a, ip_mat *b)
{
    unsigned int ch, row, col;
    ip_mat *out;

    not_null_ip_mat(a);
    not_null_ip_mat(b);

    /* faccio la verifica per vedere se sono uguali le dimensioni delle due matrici a 3 dimensioni */
    equal_dimension(a, b);
    /* siccome ho in input 2 matrici 3D che hanno la stessa dimensione , faccendo la sottrazione ottentiamo una nuova matrice 3D con le stesse dimensioni*/
    out = ip_mat_create(a->h, a->w, a->k, 0.0); /*creo la nuova matrice 3D*/
    for (ch = 0; ch < a->k; ch++)
    {
        for (row = 0; row < a->h; row++)
        {
            for (col = 0; col < a->w; col++)
                set_val(out, row, col, ch, (get_val(a, row, col, ch) - get_val(b, row, col, ch))); /* effetuo la sottrazione */
        }
    }
    compute_stats(out); /*modifico le statistiche per ogni canale della nuova matrice 3D*/
    return out;
}

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_mul_scalar(ip_mat *a, float c)
{

    unsigned int ch, row, col;
    ip_mat *out;

    not_null_ip_mat(a);

    out = ip_mat_create(a->h, a->w, a->k, 0.0); /*creo la nuova matrice 3D, che inizialmente ha tutti i valori 0.0*/
    for (ch = 0; ch < a->k; ch++)
    {
        for (row = 0; row < a->h; row++)
        {
            for (col = 0; col < a->w; col++)
                set_val(out, row, col, ch, (get_val(a, row, col, ch) * c));
        }
    }
    compute_stats(out); /*modifico le statistiche per ogni canale*/
    return out;
}

/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output.
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_add_scalar(ip_mat *a, float c)
{
    unsigned int ch, row, col;
    ip_mat *out;

    not_null_ip_mat(a);

    out = ip_mat_create(a->h, a->w, a->k, 0.0); /*creo la nuova matrice 3D*/
    for (ch = 0; ch < a->k; ch++)
    {
        for (row = 0; row < a->h; row++)
        {
            for (col = 0; col < a->w; col++)
                set_val(out, row, col, ch, (get_val(a, row, col, ch) + c));
        }
    }
    compute_stats(out); /*modifico le statistiche per ogni canale*/
    return out;
}

/* Calcola la media di due ip_mat a e b. La media si calcola per coppie delle due matrici aventi gli stessi indici
 * C[i][j][k]= (A[i][j][k]+B[i]j[k])/2
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_mean(ip_mat *a, ip_mat *b)
{
    unsigned int ch, row, col;
    ip_mat *out;

    not_null_ip_mat(a);
    not_null_ip_mat(b);
    /*faccio la verifica per vedere se sono uguali le dimensioni delle ip_mat *a e ip_mat *b */

    equal_dimension(a, b);
    /* siccome ho in input 2 matrici 3D che hanno la stessa dimensione , faccendo la media delle 2 , otteniamo una nuova matrice a 3 dimensione che ha le stesse dimensioni*/
    out = ip_mat_create(a->h, a->w, a->k, 0.0); /*creo la nuova matrice 3D*/
    for (ch = 0; ch < a->k; ch++)
    {
        for (row = 0; row < a->h; row++)
        {
            for (col = 0; col < a->w; col++)
                set_val(out, row, col, ch, (get_val(a, row, col, ch) + get_val(b, row, col, ch)) / 2);
        }
    }
    compute_stats(out); /*modifico le statistiche per ogni canale della nuova matrice a tre dimensioni*/
    return out;
}

/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/
/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 * */

ip_mat *ip_mat_to_gray_scale(ip_mat *in)
{
    unsigned int ch, row, col;
    float sup;
    ip_mat *nuova;

    not_null_ip_mat(in);

    nuova = ip_mat_create(in->h, in->w, in->k, 0.0);

    not_null_ip_mat(nuova);

    for (ch = 0; ch < nuova->k; ch++)
    {
        for (row = 0; row < nuova->h; row++)
        {
            for (col = 0; col < nuova->w; col++)
            {
                sup = mean_pixel_channel(in, row, col);
                set_val(nuova, row, col, ch, sup);
            }
        }
    }

    compute_stats(nuova);

    return nuova;
}

/* Effettua la fusione (combinazione convessa) di due immagini.
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_blend(ip_mat *a, ip_mat *b, float alpha)
{
    not_null_ip_mat(a);
    not_null_ip_mat(b);
    equal_dimension(a, b);

    if (alpha >= 0.0 && alpha <= 1.0)
    {
        ip_mat *to_blend_a, *to_blend_b, *final;
        to_blend_a = ip_mat_mul_scalar(a, alpha);
        to_blend_b = ip_mat_mul_scalar(b, 1 - alpha);
        final = ip_mat_sum(to_blend_a, to_blend_b);
        ip_mat_free(to_blend_a);
        ip_mat_free(to_blend_b);
        return final;
    }
    else
    {
        printf("Invalid alpha value\n");
        exit(1);
    }
}

/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_brighten(ip_mat *a, float bright)
{
    not_null_ip_mat(a);

    return ip_mat_add_scalar(a, bright);
}

/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 *
 * out = a + gauss_noise*amount
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 * */
ip_mat *ip_mat_corrupt(ip_mat *a, float amount)
{
    not_null_ip_mat(a);

    if (amount >= MIN_PIXEL_FLOAT && amount <= MAX_PIXEL_FLOAT)
    {
        unsigned int ch, row, col;
        ip_mat *out;
        ip_mat *temp = ip_mat_create(a->h, a->w, a->k, 0.0);
        float sup = 0.0;

        for (ch = 0; ch < a->k; ch++)
        {
            for (row = 0; row < a->h; row++)
            {
                for (col = 0; col < a->w; col++)
                {
                    sup = get_normal_random(0.0, amount / 2);
                    set_val(temp, row, col, ch, sup);
                }
            }
        }
        out = ip_mat_sum(a, temp);
        ip_mat_free(temp);

        return out;
    }
    else
    {
        printf("Provide a valid amount value!\n");
        exit(1);
    }
}

/**** PARTE 3: CONVOLUZIONE E FILTRI *****/

/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a".
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat *ip_mat_convolve(ip_mat *a, ip_mat *f)
{
    ip_mat *padded, *out;
    channel_t fil_ch;
    unsigned int ch, row, col, pad_h, pad_w;
    not_null_ip_mat(a);
    not_null_ip_mat(f);
    pad_h = (f->h - 1) / 2;
    pad_w = (f->w - 1) / 2;
    /* inizializza matrici per il calcolo della convolve */
    padded = ip_mat_padding(a, pad_h, pad_w);
    out = ip_mat_create(a->h, a->w, a->k, 0.0);

    for (ch = 0; ch < a->k; ch++)
    {
        /* questa operazione assicura che vengano applicati i primi a->k canali del filtro all'immagine, e se il filtro non ha canali sufficienti si applica sempre l'ultimo canale del filtro disponibile */
        if (ch < f->k)
            fil_ch = get_channel(f, ch);
        for (row = 0; row < a->h; row++)
            for (col = 0; col < a->w; col++)
                set_val(out, row, col, ch, convolve_channel(get_channel(padded, ch), fil_ch, row, col));
    }
    /* libera la matrice temporanea */
    ip_mat_free(padded);
    return out;
}

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
ip_mat *ip_mat_padding(ip_mat *a, unsigned int pad_h, unsigned int pad_w)
{
    ip_mat *out;

    not_null_ip_mat(a);

    out = ip_mat_create(a->h + 2 * pad_h, a->w + 2 * pad_w, a->k, 0.0);
    ip_mat_puts(out, a, pad_h, pad_w, COMPUTE_STATS);

    return out;
}

/* Crea un filtro di sharpening */
ip_mat *create_sharpen_filter()
{
    /* crea un filtro di sharpening semplice, a un canale */
    ip_mat *out = ip_mat_create(3, 3, 1, 0.0);

    set_val(out, 0, 1, 0, -1.0);
    set_val(out, 1, 0, 0, -1.0);
    set_val(out, 1, 1, 0, 5.0);
    set_val(out, 1, 2, 0, -1.0);
    set_val(out, 2, 1, 0, -1.0);
    return out;
}

/* Crea un filtro per rilevare i bordi */
ip_mat *create_edge_filter()
{
    /* crea un filtro di edge semplice, a un canale */
    ip_mat *out = ip_mat_create(3, 3, 1, -1.0);

    set_val(out, 1, 1, 0, 8.0);
    return out;
}

/* Crea un filtro per aggiungere profondità */
ip_mat *create_emboss_filter()
{
    /* crea un filtro di emboss semplice, a un canale */
    ip_mat *out = ip_mat_create(3, 3, 1, 1.0);

    set_val(out, 0, 0, 0, -2.0);
    set_val(out, 0, 1, 0, -1.0);
    set_val(out, 0, 2, 0, 0.0);
    set_val(out, 1, 0, 0, -1.0);
    set_val(out, 2, 0, 0, 0.0);
    set_val(out, 2, 2, 0, 2.0);
    return out;
}

/* Crea un filtro medio per la rimozione del rumore */
ip_mat *create_average_filter(unsigned int h, unsigned int w, unsigned int k)
{
    float c = 1.0 / (w * h);
    /* crea un filtro average a k canali */
    ip_mat *out = ip_mat_create(h, w, k, c);

    return out;
}

/* Crea un filtro gaussiano per la rimozione del rumore */
ip_mat *create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma)
{
    if (h >= 3 && w >= 3 && k > 0 && sigma > 0)
    {
        unsigned int cx, cy, row, col, ch;
        float *sums = (float *)malloc(sizeof(float) * k);
        ip_mat *out;
        out = ip_mat_create(h, w, k, 0.0);

        cx = w / 2;
        cy = h / 2;
        for (ch = 0; ch < k; ch++)
        {
            sums[ch] = 0.0;
            for (row = 0; row < h; row++)
            {
                for (col = 0; col < w; col++)
                {
                    int x, y;
                    float val;
                    x = row - cx;
                    y = col - cy;
                    val = (1. / (2. * PI * sigma * sigma)) * exp(-(x * x + y * y) / (2. * sigma * sigma)); /*il valore da inserire , calcolato secondo la formula data*/
                    set_val(out, row, col, ch, val);                                                       /*inserisco il valore nella sua posizione*/
                    sums[ch] += get_val(out, row, col, ch);
                }
            }
        }

        /* normalizzazione */
        for (ch = 0; ch < k; ch++)
        {
            for (row = 0; row < h; row++)
            {
                for (col = 0; col < w; col++)
                {
                    float new_val;
                    new_val = get_val(out, row, col, ch) / sums[ch];
                    set_val(out, row, col, ch, new_val);
                }
            }
        }
        free(sums);
        return out;
    }
    else
    {
        printf("Passed parameters are invalid , please enter other parameters");
        exit(1);
    }
}

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
void rescale(ip_mat *t, float new_max)
{
    unsigned int ch, row, col;
    not_null_ip_mat(t);
    for (ch = 0; ch < t->k; ch++)
    {
        for (row = 0; row < t->h; row++)
        {
            for (col = 0; col < t->w; col++)
            {
                float val, max, min;
                max = t->stat[ch].max;
                min = t->stat[ch].min;
                val = (get_val(t, row, col, ch) - min) / (max - min) * new_max;
                set_val(t, row, col, ch, val);
            }
        }
    }
}

/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.
 *
 * Il risultato dell'operazione si salva in t
 * */
void clamp(ip_mat *t, float low, float high)
{
    unsigned int ch, row, col;
    not_null_ip_mat(t);

    for (ch = 0; ch < t->k; ch++)
    {
        for (row = 0; row < t->h; row++)
        {
            for (col = 0; col < t->w; col++)
                set_val(t, row, col, ch, restrict_val(get_val(t, row, col, ch), low, high));
        }
    }
}

/**** METODI GIA' IMPLEMENTATI ****/
/* Genera dei numeri casuali con distribuzione Normale
 * https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
 * caratterizzata da media = mu e deviazione  = std
 * */
float get_normal_random(float media, float std)
{

    float y1 = ((float)(rand()) + 1.) / ((float)(RAND_MAX) + 1.);
    float y2 = ((float)(rand()) + 1.) / ((float)(RAND_MAX) + 1.);
    float num = cos(2 * PI * y2) * sqrt(-2. * log(y1));

    return media + num * std;
}

/* Converte una Bitmap in una ip_mat*/
ip_mat *bitmap_to_ip_mat(Bitmap *img)
{
    if (img)
    {
        unsigned int i = 0, j = 0;

        unsigned char R, G, B;

        unsigned int h = img->h;
        unsigned int w = img->w;

        ip_mat *out = ip_mat_create(h, w, 3, 0);

        for (i = 0; i < h; i++) /* rows */
        {
            for (j = 0; j < w; j++) /* columns */
            {
                bm_get_pixel(img, j, i, &R, &G, &B);
                set_val(out, i, j, 0, (float)R);
                set_val(out, i, j, 1, (float)G);
                set_val(out, i, j, 2, (float)B);
            }
        }

        compute_stats(out);

        return out;
    }
    else
    {
        printf("Insert a valid image\n");
        exit(1);
    }
}

/* Converte una ip_mat in una bitmap*/
Bitmap *ip_mat_to_bitmap(ip_mat *t)
{
    not_null_ip_mat(t);
    {
        Bitmap *b = bm_create(t->w, t->h);

        unsigned int i, j;
        for (i = 0; i < t->h; i++)
        {
            for (j = 0; j < t->w; j++)
            {
                bm_set_pixel(b, j, i, (unsigned char)restrict_val(get_val(t, i, j, 0), MIN_PIXEL_FLOAT, MAX_PIXEL_FLOAT),
                             (unsigned char)restrict_val(get_val(t, i, j, 1), MIN_PIXEL_FLOAT, MAX_PIXEL_FLOAT),
                             (unsigned char)restrict_val(get_val(t, i, j, 2), MIN_PIXEL_FLOAT, MAX_PIXEL_FLOAT));
            }
        }
        return b;
    }
}

/* Visualizza i dati stampando in ordine le matrici rispetto
 * la terza dimensione.
 * Prima stamperemo t->data[...][...][0] poi t->data[...][...][1] ...
 * */
void ip_mat_show(ip_mat *t)
{
    unsigned int i, l, j;
    not_null_ip_mat(t);
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
    not_null_ip_mat(t);

    compute_stats(t);

    for (k = 0; k < t->k; k++)
    {
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}
