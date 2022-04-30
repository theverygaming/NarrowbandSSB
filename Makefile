CFLAGS = -Wall -Wshadow -O3 -g
LDLIBS = -lsndfile -lm -lfftw3 -lfftw3f -lvolk
CC = g++

all: main

main: modulate.o demodulate.o
	$(CC) modulate.o $(LDLIBS) -o build/modulate
	$(CC) demodulate.o $(LDLIBS) -o build/demodulate

clean:
	$(RM) *.o
	$(RM) build/modulate
	$(RM) build/demodulate
