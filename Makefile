CC = cc
CFLAGS = -Wall -O2 -Ilib -Ilib/laspack -Ilib/xc
LD = cc
LFLAGS = -lm 
LASPACKOFILES = ./lib/laspack/eigenval.o \
								./lib/laspack/errhandl.o \
								./lib/laspack/factor.o \
								./lib/laspack/itersolv.o \
								./lib/laspack/matrix.o \
								./lib/laspack/mlsolv.o \
								./lib/laspack/operats.o \
								./lib/laspack/precond.o \
								./lib/laspack/qmatrix.o \
								./lib/laspack/rtc.o \
								./lib/laspack/vector.o

OFILES = ./src/calculation.o \
         ./src/initialize.o \
         ./src/main.o \
         ./src/norm.o \
         ./src/construction.o

all: $(OFILES) $(LASPACKOFILES)
	$(LD) $(LFLAGS) -o gas $(OFILES) $(LASPACKOFILES)

%.o: %.c
	$(CC) $(CFLAGS) -o $*.o -c $*.c
	
clean:
	rm -f $(OFILES) $(LASPACKOFILES)
