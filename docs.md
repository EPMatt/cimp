# iplib - Documentazione di progetto

##### Gruppo n.20 - MVP

##### v1.0.0 (31/05/2020)

## Indice

- [Introduzione](#introduzione)
- [Note importanti](#note-importanti)
- [Funzioni modificate](#funzioni-modificate)
- [Funzioni di supporto](#funzioni-di-supporto)
  - [Generiche](#generiche)
  - [Operazioni su ip_mat](#operazioni-su-ip_mat)
  - [Operazioni su channel_t](#operazioni-su-channel_t)
- [Makefile](#makefile)
- [Testing](#testing)
- [Gestione degli errori](#gestione-degli-errori)
- [Membri del gruppo](#membri-del-gruppo)

## Introduzione

Questo documento raggruppa alcune note importanti in merito allo sviluppo della libreria iplib. Il progetto sarà disponibile su GitHub a [questo link](https://github.com/EPMatt/c-2020).

## Note importanti

L'implementazione del campo `data` della struttura `ip_mat` prevede un array tridimensionale ad indici frastagliati del tipo *canale x righe x colonne*. In questo modo si ha una netta separazione tra i diversi canali della matrice, semplificando numerose operazioni strettamente operanti su canali.

A tale scopo si è definito il tipo `channel_t`, corrispondente ad un singolo canale di una `ip_mat`.

```c
    typedef struct
    {
        float **data;
        unsigned int h;
        unsigned int w;
    } channel_t;
```

Strutture di tipo `channel_t` vengono utilizzate nella libreria per fornire un accesso diretto ad un singolo canale di una matrice di pixel.

## Funzioni modificate

Tra le funzioni preimplementate alcune sono state modificate come indicato in seguito.

*`get_val` / `set_val`: modificato l’ordine delle dimensioni all'accesso dei dati nel campo `data` della struttura `ipmat` passata come parametro.
*`ip_mat_to_bitmap` / `bitmap_to_ip_mat` e `ip_mat_show`/`ip_mat_show_stats`:  aggiunto controllo se il puntatore a `ip_mat` o a `bitmap` non sia `NULL`.
*`ip_mat_to_bitmap` nella conversione dell'immagine i valori dei pixel recuperati dalla `ip_mat` passata come parametro vengono ristretti nell'intervallo di valori consentito per un pixel. Per svolgere tale compito si utilizza la funzione ausiliaria `restrict_val`.

## Funzioni di supporto

Per la realizzazione della libreria si è deciso di implementare alcune funzioni ausiliarie per migliorare la strutturazione del codice (Maggiori dettagli delle funzioni in ip_lib.c).
Tali funzioni sono qui di seguito elencate.

### Generiche

#### `float restrict_val(float val, float low, float high);`

Restringe il valore fornito all'interno del range low...high, estremi inclusi.

### Operazioni su ip_mat

#### `float min(ip_mat *a, int k);`

Funzione ausiliaria min che trova il minimo di un determinato canale k, riceve in input una ip_mat a e un canale k e restituisce un float che è il minimo valore di pixel trovato nel canale k.

#### `float max(ip_mat *a, int k);`

Funzione ausiliaria che calcola il valore massimo di un determinato canale k. Funzionamento analogo al sottoprogramma precedente.

#### `float mean(ip_mat *a, int k);`

Funzione ausiliaria che calcola la media degli elementi che si trovano in un determinato canale k. Funzionamento analogo al sottoprogramma precedente.

#### `float mean_pixel_channel(ip_mat *a, unsigned int i, unsigned int j);`

Calcola la media per il pixel (i,j) sui k canali della matrice passata come parametro.

#### `void ip_mat_puts(ip_mat *dest, ip_mat const *source, unsigned int row, unsigned int col, int do_compute_stats);`

Posiziona il contenuto della matrice `source` nella matrice `dest`, partendo dalla posizione specificata.
Per tutti i canali di `source` (e massimo il numero di canali di `dest`) inizia posizionando la cella 0,0 del canale di `source` alla posizione `row`,`col` del canale di `dest`.
Quindi copia il contenuto di `source` fino a raggiungere il limite del canale di `dest`, e soltanto se `row`,`col` sono indici validi per la matrice `dest`.
Opzionalmente, se il parametro `do_compute_stats` è pari a `COMPUTE_STATS`, ricalcola le statistiche sulla matrice fornita.

#### `void not_null_ip_mat(ip_mat *a);`

Restituisce un messaggio di errore nel caso in cui la matrice passata sia `NULL`.

#### `void equal_dimension(ip_mat *a, ip_mat *b);`

Controlla se due `ip_mat` hanno le stesse dimensioni oppure no.
Esce dal programma nel momento in cui una delle dimensioni é diversa.

### Operazioni su channel_t

#### `channel_t get_channel(ip_mat const *a, unsigned int ch);`

Data una `ip_mat`, ottieni una struttura canale, con riferimento al ch-esimo canale della matrice. Si tratta appositamente di un riferimento, non di una copia del canale della matrice fornita come parametro.

#### `void channel_puts(channel_t const dest, channel_t const source, unsigned int row, unsigned int col);`

Posiziona il contenuto del canale `source` nel canale `dest`. Inizia posizionando la cella 0,0 del canale `source` alla posizione `row`,`col` del canale `dest`.
Quindi copia il contenuto di `source` fino a raggiungere il limite del canale `dest`, e soltanto se `row`,`col` sono indici validi per il canale `dest`.

#### `float convolve_channel(channel_t ch, channel_t filter, unsigned int start_h, unsigned int start_w);`

Calcola la somma di prodotti tra il kernel fornito e il canale fornito, partendo dalla posizione (start_h,start_w) del canale.

## Makefile

Il Makefile fornito con questo progetto fornisce diverse funzionalità aggiuntive rispetto alla semplice compilazione per l'ambiente di production.
Si forniscono le seguenti recipe, richiamabili con il comando `make`:
*`build` (default) compila i sorgenti per l'ambiente di production;
*`test` compila i sorgenti per l'ambiente di testing;
*`test-run` avvia il programma di test attraverso l'utility CLI `valgrind`;
*`clean` elimina tutti i file risultanti dalla compilazione.

Make è in grado di determinare automaticamente se i file compilati sono stati generati per testing o production, e all'occorrenza effettuare un clean prima di ricompilare i sorgenti. Questo meccanismo è possibile grazie all'utilizzo del file esterno `.lastmake`, compilato dalle recipe `build` e `test`.

Allo stesso modo, non è possibile eseguire la recipe `test-run` se i sorgenti sono stati compilati per un ambiente "production". Sarà necessario prima eseguire la recipe `test` e successivamente `test-run`.

## Testing

Si ritiene opportuno includere nella consegna `test_iplib.c`, il file utilizzato durante lo sviluppo per effettuare test di funzionamento della libreria.
Tale programma esegue un numero prestabilito di round di test.
Ciascun round è caratterizzato dalla generazione di una `ipmat` di dimensione casuale (limitati superiormente) e l'esecuzione di tutte le funzioni della libreria in successione. I risultati del round vengono memorizzati in un array, che in caso di necessità sarà facilmente ispezionabile. Il programma è costruito per essere eseguito con l'utility CLI valgrind, e verificare l'assenza di memory leaks o altre operazioni errate sulla memoria, per ciascuna funzione della libreria. Non rappresenta un flusso normale di esecuzione e non è in alcun modo correlato con il programma `main_iplib.c`.

## Gestione degli errori

L'approccio adottato per la gestione degli errori si basa sul principio di limitare il più possibile la propagazione di un errore dovuto, ad esempio al passaggio ad un sottoprogramma di parametri non adeguati o all'impossibilità di allocare un'area di memoria necessaria per lo svolgimento del task richiesto. Per questo motivo ciascuna funzione si assicura, prima si svolgere qualsiasi operazione, di aver ricevuto parametri corretti. In caso negativo, l'esecuzione del programma sarà interrotta e un messaggio di errore sarà mostrato all'utente.

## Membri del gruppo

-[Ina Popescu](https://github.com/Ina-pps) (matricola ***REMOVED***)
-[Matteo Agnoletto](https://github.com/EPMatt) (matricola ***REMOVED***)
-[Lorenzo Armando Donatelli](https://github.com/Donnyz) (matricola: ***REMOVED***)
