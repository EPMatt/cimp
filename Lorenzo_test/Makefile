ip_lib.o: ip_lib.c ip_lib.h
	gcc ip_lib.c -o ip_lib.o --ansi --pedantic -Wall -c -lm

test_bmp: test_bmp.c bmp.o bmp.h
	gcc test_bmp.c bmp.o -o test_bmp -Wall -lm

bmp.o: bmp.c bmp.h
	gcc bmp.c -o bmp.o -Wall -lm -c 

c:
	clear
