CC=gcc
CFLAGS= -Ofast 

calc: main.c bessel_mod.c 
	$(CC) -o calc bessel_mod.c -fopenmp main.c -lm -lgsl -lgslcblas
