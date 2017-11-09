# Compiler flags
#  -g		adds debugging information
#  -Wall	turns on most compiler warnings
#  -fbounds-check
CFLAGS = -g -Wall

HEADERS = percolation.h
OBJECTS = percolation.o linked_list.o

# Build target executable
default: main

main.o: main.c $(OBJECTS) $(HEADERS)
	gcc -c main.c -o main.o

percolation.o: percolation.c $(HEADERS)
	gcc -c percolation.c -o percolation.o

linked_list.o: linked_list.c $(HEADERS)
	gcc -c linked_list.c -o linked_list.o

main: main.o
	gcc -o main main.o percolation.o linked_list.o
	
clean:
	-rm -f *.o
	-rm -f main
