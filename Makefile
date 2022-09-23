CC=g++

all:
	$(CC) -std=c++11 liddriven.cpp -o liddriven.o

clean:
	rm -r liddriven.o 
