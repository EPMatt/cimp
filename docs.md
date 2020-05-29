# iplib - Documentazione di progetto

- [iplib - Documentazione di progetto](#iplib---documentazione-di-progetto)
  - [Introduzione](#introduzione)
  - [Note importanti](#note-importanti)
  - [Funzioni modificate](#funzioni-modificate)
  - [Funzioni di supporto](#funzioni-di-supporto)
  - [Makefile](#makefile)

## Introduzione
## Note importanti
Per la realizzazione del file ip_lib.c si è preferito, come visualizzazione delle matrici ip_mat l’ordine delle dimensioni k h w, cioè i campi che puntano ad array  multidimensionali frastagliati di dimensione h x w.
## Funzioni modificate
Tra le funzioni che sono state fornite sono state modificate le seguenti nel modo indicato sotto: 
* get_val / set_val:Modificato l’ordine delle dimensioni come detto nelle note importanti
* ip_mat_to_bitmap / bitmat_to_ip_mat e ip_mat_show/ip_mat_show_stats:  Aggiunto controllo se il puntatore a ip_mat o a bitmap non punti a null


## Funzioni di supporto
Per la realizzazione della libreria si è preso la decisione di implementare e arricchire il file ip_lib.c con delle funzioni ausiliarie in modo da rendere il codice più leggibile e ridurre la riscrittura del medesimo codice(stasera non mi viene il termine).
(Maggiori dettagli delle funzioni in ip_lib.c)

Tali funzioni sono qui di seguito elencate e suddivise in tipo restituito e scopo:
    • Float
        ◦ min/max/mean (ip_mat* t, unsigned int k):
              Funzioni matematiche che data una ip_mat e un unsigned int che rappresenta il campo calcolano rispettivamente minimo, massimo e media di un campo.

        ◦ restrict_val (float val, float low, float high)
              Funzione che dato un valore val, controlla se esso è compreso tra un valore minimo (low) e un valore massimo (high). 
              Se tale valore supera questo range, verrà sostituito con il valore dell’estremo che esso supera.
              
        ◦ mean_pixel_channel(ip_mat* t, unsigned int i, unsigned int j):
		Funzione che data la posizione di un pixel indicata con i valori i e j, calcola la media
		dei pixel in posizione i e j dei tre canali.
      
## Makefile