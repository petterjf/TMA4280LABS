CC=gcc
CFLAGS=-lm -Wall

all: num_pi

num_pi: num_pi.c zeta.c
	$(CC) $(CFLAGS) -o num_pi.out num_pi.c zeta.c

utest: utest.c zeta.c
	$(CC) $(CFLAGS) -o utest.out utest.c zeta.c

vtest: vtest.c zeta.c
	$(CC) $(CFLAGS) -o vtest.out vtest.c zeta.c