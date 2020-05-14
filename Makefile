bmp.o:bmp.c bmp.h
		gcc bmp.c -c -Wall -lm

test_bmp:test_bmp.c bmp.o bmp.h
		gcc test_bmp.c bmp.o -o test_bmp -Wall -lm

clean:
		rm -f *.o test_bmp
