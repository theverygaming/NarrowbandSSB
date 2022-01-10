CFLAGS = -Wall -Wshadow -O3 -g
LDLIBS = -lsndfile -lm -lfftw3 -lvolk
CC = gcc

main: decimate.o filter.o
	$(CC) decimate.o filter.o $(LDLIBS) -o build/decimate

clean:
	$(RM) *.o
	$(RM) build/decimate
