OBJECTS = main_iplib.o bmp.o ip_lib.o

BUILD_FLAGS= -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
TEST_FLAGS= -Wall -lm -g -O1
MAKE_ENV=build
-include .lastmake

build:
		@echo "Building sources..."
ifneq ($(LAST_MAKE),build)
		@echo "Ooops! Last make was not a build one, must delete all generated files first!"
		@make clean
		@echo "Ready to build sources"
endif
		@echo 'LAST_MAKE=build' > .lastmake
		@make -e main_ip_lib
		@echo "...done!"

test:
		@echo "Building sources for testing..."
ifneq ($(LAST_MAKE),test)
		@echo "Ooops! Last make was not a test one, must delete all generated files first!"
		@make clean
		@echo "Ready to build sources for testing"
endif
		@echo 'LAST_MAKE=test' > .lastmake
		@make -e MAKE_ENV=test main_ip_lib test_iplib
		@echo "...done!"

test-run:
ifneq ($(LAST_MAKE),test)
		@echo "Ooops! Last make was not a test one, you can't run memory tests on a production build."
		@echo "Run make test first, then run this command again."
else
		valgrind -v --leak-check=full --track-origins=yes ./test_iplib 
endif
	

clean:
		@echo "Cleaning compiled files..."
		@rm -f test_bmp main_iplib test_iplib $(OBJECTS) .lastmake
		@echo "...done!"

test_iplib: test_iplib.c bmp.o ip_lib.o ip_lib.h bmp.h
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

.PHONY: build test test-run clean