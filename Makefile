main_ip_lib: main_iplib.o bmp.o ip_lib.o
		gcc main_iplib.o ip_lib.o bmp.o -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra -o main_iplib
		
main_iplib.o: main_iplib.c ip_lib.h bmp.h
		gcc main_iplib.c -c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

ip_lib.o: ip_lib.c ip_lib.h bmp.h
		gcc ip_lib.c -c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

bmp.o:bmp.c bmp.h
		gcc bmp.c -c -Wall -lm

test_bmp:test_bmp.c bmp.o bmp.h
		gcc test_bmp.c bmp.o -o test_bmp -Wall -lm

clean:
		rm -f *.o test_bmp main_iplib