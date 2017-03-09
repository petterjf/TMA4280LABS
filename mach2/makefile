CC=gcc
CFLAGS=-Wall -lm -fopenmp

all: num_pi

num_pi: num_pi.c mach.c
	$(CC) $(CFLAGS) -o num_pi.out num_pi.c mach.c