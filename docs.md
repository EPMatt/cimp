# iplib - Documentazione di progetto gruppo MVP 

- [iplib - Documentazione di progetto](#iplib---documentazione-di-progetto)
  - [Introduzione](#introduzione)
  - [Note importanti](#note-importanti)
  - [Funzioni modificate](#funzioni-modificate)
  - [Funzioni di supporto](#funzioni-di-supporto)
  - [Makefile](#makefile)

## Introduzione
Spiegazioni basilari funzioni ausiliarie gruppo MVP
## Note importanti
Per la realizzazione del file ip_lib.c si è preferito, come visualizzazione delle matrici ip_mat l’ordine delle dimensioni k h w, cioè i campi che puntano ad array  multidimensionali frastagliati di dimensione h x w.

Creazione della struttura/tipo channel_t corrispondente ad una matrice di float, utilizzata per alcune funzioni:

	typedef struct
	{
		float **data;
		unsigned int h;
		unsigned int w;
	} channel_t;

## Funzioni modificate
Tra le funzioni che sono state fornite sono state modificate le seguenti nel modo indicato sotto: 

- get_val / set_val:Modificato l’ordine delle dimensioni come detto precedentemente
- ip_mat_to_bitmap / bitmat_to_ip_mat e ip_mat_show/ip_mat_show_stats:  Aggiunto controllo se il puntatore a ip_mat o a bitmap non punti a null

- ip_mat_to_bitmap
	prima di convertire l'immagine effetuiamo un controllo dei valori con la funzione restrict_val


## Funzioni di supporto
Per la realizzazione della libreria si è preso la decisione di implementare e arricchire il file ip_lib.c con delle funzioni ausiliarie in modo da rendere il codice più leggibile e ridurre la riscrittura del medesimo codice(stasera non mi viene il termine).
(Maggiori dettagli delle funzioni in ip_lib.c)

Tali funzioni sono qui di seguito elencate e suddivise in tipo restituito e scopo:

- ## Float
                                                                                                                
		min(ip_mat *a, int k);
			funzione aux min che trova il minimo di un determinato canale k, riceve in input una ip_mat a e un canale k e restituisce un float che è il minimo


		max(ip_mat *a, int k);
			funzione aux che mi calcola il valore massimo di un determinato canale k


		mean(ip_mat *a, int k);
			mean è la funzione aux che mi permette di calcolare la media degli elementi che si trovano in un determinato canale k*/


		restrict_val(float val, float low, float high);
			restringi il valore fornito all'interno del range low...high, estremi inclusi 
		

		mean_pixel_channel(ip_mat *a, unsigned int i, unsigned int j);
			 Calcola la media per ogni pixel sui tre canali, prendendo un indice colonna e un indice riga.

		convolve_channel(channel_t ch, channel_t filter, unsigned int start_h, unsigned int start_w);
			calcola la somma di prodotti tra il kernel fornito e il canale fornito, partendo dalla posizione (start_h,start_w) del canale

- ## Channel_t
		get_channel(ip_mat const *a, unsigned int ch);
			data una ip_mat, ottieni una struttura canale, con riferimento al ch-esimo canale 
 
- ## Void
		channel_puts(channel_t const dest, channel_t const source, unsigned int row, unsigned int col);
			posiziona il contenuto del canale source nel canale dest, posizionando la cella 0,0 del canale source alla posizione row,col del canale dest
			copia il contenuto di source fino a raggiungere il limite del canale dest, e soltanto se row,col sono indici validi per il canale dest
 
		ip_mat_puts(ip_mat *dest, ip_mat const *source, unsigned int row, unsigned int col, int do_compute_stats);
			posiziona il contenuto della matrice source nella matrice dest, partendo dalla posizione specificata
 			per tutti i canali di source (massimo il numero di canali di dest) posiziona la cella 0,0 del canale di source alla posizione row,col del canale di dest
 			copia il contenuto di source fino a raggiungere il limite del canale di dest, e soltanto se row,col sono indici validi per la matrice dest

		not_null_ip_mat(ip_mat *a);
			restituisce un messaggio di errore nel caso in cui la matrice passata sia NULL

		two_not_null_ip_mat(ip_mat *a, ip_mat *b);
 			restituisce un messaggio di errore nel caso in cui una o più ip_mat siano null

		equal_dimension(ip_mat *a, ip_mat *b);
			Controlla se due ip_mmat hanno le stesse dimensioni oppre no
			esce dal programma nel momento in cui una delle dimensioni é diversa





## Makefile