CFLAGS = -Wall -Wshadow -O3 -g
LDLIBS = -lsndfile -lm -lfftw3 -lvolk
CC = g++

main: modulate.o demodulate.o filter.o
	$(CC) modulate.o filter.o $(LDLIBS) -o build/modulate
	$(CC) demodulate.o filter.o $(LDLIBS) -o build/demodulate

clean:
	$(RM) *.o
	$(RM) build/modulate
	$(RM) build/demodulate
