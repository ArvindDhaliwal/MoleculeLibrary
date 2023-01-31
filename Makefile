all:test 

test: test.o libmol.so 
	clang test.o -L. -lmol -lm -o test

test.o: test.c mol.h
	clang -c -Wall -std=c99 test.c 

mol.o: mol.c
	clang -c -Wall -std=c99 -fPIC mol.c

libmol.so: mol.o
	clang mol.o -shared -o libmol.so

clean: 
	rm -f *.o mol test
