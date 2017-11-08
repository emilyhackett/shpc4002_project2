# Compiler flags
#  -g		adds debugging information
#  -Wall	turns on most compiler warnings
#  -fbounds-check
CFLAGS = -g -Wall

# Build target executable
default: main

main.o: main.c percolation.o percolation.h
	gcc -c main.c -o main.o

percolation.o: percolation.c percolation.h
	gcc -c percolation.c -o percolation.o

main: main.o
	gcc -o main main.o percolation.o
	
clean:
	-rm -f main.o
	-rm -f main
	-rm -f percolation.o
