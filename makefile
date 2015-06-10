CC = gcc
CCFLAGS = 

CC_COMPILE  =  $(CC) $(CCFLAGS) -I/usr/local/include -c -O3
CC_LOAD     = $(CC) $(CCFLAGS) -L/usr/local/lib 

.c.o:
	$(CC_COMPILE) $*.c

EXTRALIBS =  -lm -lgsl -lgslcblas

EXE = iharmony
all: $(EXE)

SRCS = bessel_mod.c main.c integrate.c
OBJS = bessel_mod.o main.o integrate.o

$(OBJS) : $(INCS) makefile


$(EXE): $(OBJS) $(INCS) makefile
	$(CC_LOAD) $(OBJS) $(EXTRALIBS) -o $(EXE)

clean:
	/bin/rm *.o
	/bin/rm $(EXE)