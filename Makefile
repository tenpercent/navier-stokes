CC = cc
CFLAGS = -Wall -O2 -Ilib -Ilib/laspack
LD = cc
LFLAGS = -lm 
LASPACKOFILES = ./lib/laspack/qmatrix.o \
                ./lib/laspack/vector.o \
                ./lib/laspack/itersolv.o \
                ./lib/laspack/rtc.o \
                ./lib/laspack/matrix.o \
                ./lib/laspack/errhandl.o \
                ./lib/laspack/operats.o \
                ./lib/laspack/eigenval.o
OFILES = ./src/calculation.o \
         ./src/initialize.o \
         ./src/main.o

all: $(OFILES) $(LASPACKOFILES)
	$(LD) $(LFLAGS) -o gas $(OFILES) $(LASPACKOFILES)
%.o: %.c
	$(CC) $(CFLAGS) -o $*.o -c $*.c
clean:
	rm -f $(OFILES) $(LASPACKOFILES)

