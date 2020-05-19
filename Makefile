objects = main_iplib.o bmp.o ip_lib.o

build_flags= -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

build: main_ip_lib

clean:
		rm -f test_bmp main_iplib $(objects)

main_ip_lib: $(objects) ip_lib.h bmp.h 
		gcc $(objects) $(build_flags) -o main_iplib

main_iplib.o: main_iplib.c ip_lib.h bmp.h
		gcc main_iplib.c -c $(build_flags)

ip_lib.o: ip_lib.c ip_lib.h bmp.h
		gcc ip_lib.c -c $(build_flags)

bmp.o:bmp.c bmp.h
		gcc bmp.c -c -Wall -lm

test_bmp:test_bmp.c bmp.o bmp.h
		gcc test_bmp.c bmp.o -o test_bmp -Wall -lm

.PHONY: clean build