OBJECTS = main_iplib.o bmp.o ip_lib.o

BUILD_FLAGS= -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
TEST_FLAGS= -Wall -lm -g -O1
MAKE_ENV=build

build:
		@echo "\nBuilding sources...\n"
		@make -e MAKE_ENV=build main_ip_lib
		@echo "\n...done!\n"

test:
		@echo "\nBuilding sources for testing...\n"
		@make -e MAKE_ENV=test main_ip_lib test_ip_lib
		@echo "\n...done!\n"

clean:
		@echo "\nCleaning compiled files...\n"
		@rm -f test_bmp main_iplib $(OBJECTS)
		@echo "\n...done!\n"

test_ip_lib: bmp.o ip_lib.o ip_lib.h bmp.h
	gcc test_iplib.c bmp.o ip_lib.o $(TEST_FLAGS) -o test_iplib

main_ip_lib: $(OBJECTS) ip_lib.h bmp.h
ifeq ($(MAKE_ENV),test)
		gcc $(OBJECTS) $(TEST_FLAGS) -o main_iplib
else
		gcc $(OBJECTS) $(BUILD_FLAGS) -o main_iplib
endif

main_iplib.o: main_iplib.c ip_lib.h bmp.h
ifeq ($(MAKE_ENV),test)
		gcc main_iplib.c -c $(TEST_FLAGS)
else
		gcc main_iplib.c -c $(BUILD_FLAGS)
endif

ip_lib.o: ip_lib.c ip_lib.h bmp.h
ifeq ($(MAKE_ENV),test)
		gcc ip_lib.c -c $(TEST_FLAGS)
else
		gcc ip_lib.c -c $(BUILD_FLAGS)
endif

bmp.o:bmp.c bmp.h
		gcc bmp.c -c -Wall -lm

test_bmp:test_bmp.c bmp.o bmp.h
		gcc test_bmp.c bmp.o -o test_bmp -Wall -lm

.PHONY: clean build test