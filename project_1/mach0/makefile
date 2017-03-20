CC=gcc
CFLAGS=-lm -Wall

all: num_pi

num_pi: num_pi.c mach.c
	$(CC) $(CFLAGS) -o num_pi.out num_pi.c mach.c

utest: utest.c mach.c
	$(CC) $(CFLAGS) -o utest.out utest.c mach.c

vtest: vtest.c mach.c
	$(CC) $(CFLAGS) -o vtest.out vtest.c mach.c