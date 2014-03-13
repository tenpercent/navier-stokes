BUILDDIR = build
CC = cc
CFLAGS = -Wall -O2 -Ilib -Ilib/laspack -Ilib/xc
LD = cc
LFLAGS = -lm
LASPACKONAMES = eigenval.o \
                errhandl.o \
                factor.o \
                itersolv.o \
                matrix.o \
                operats.o \
                precond.o \
                qmatrix.o \
                rtc.o \
                vector.o
LASPACKOFILES = $(foreach fname,$(LASPACKONAMES),$(BUILDDIR)/laspack/$(fname))
ONAMES = calculation.o \
         construction.o \
         initialize.o \
         main.o \
         norm.o
OFILES = $(foreach fname,$(ONAMES),$(BUILDDIR)/program/$(fname))

%/create-stamp:
	mkdir -p $*
	touch $@

gas: $(OFILES) $(LASPACKOFILES)
	$(LD) $(LFLAGS) -o gas $(OFILES) $(LASPACKOFILES)

$(BUILDDIR)/laspack/%.o: $(BUILDDIR)/laspack/create-stamp lib/laspack/%.c
	$(CC) $(CFLAGS) -o $@ -c lib/laspack/$*.c

$(BUILDDIR)/program/%.o: $(BUILDDIR)/program/create-stamp src/%.c
	$(CC) $(CFLAGS) -o $@ -c src/$*.c

clean:
	rm -rf $(BUILDDIR)

distclean: clean
	rm -f gas
