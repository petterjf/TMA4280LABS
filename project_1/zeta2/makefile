CC=gcc
CFLAGS=-Wall -fopenmp -lm

all: num_pi

num_pi: num_pi.c zeta.c
	$(CC) $(CFLAGS) -o num_pi.out num_pi.c zeta.c