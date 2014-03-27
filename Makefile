BUILDDIR = build
DEFINES =
CC = cc
CFLAGS = -Wall -g -O2 -Ilib -Ilib/laspack
CPPFLAGS = $(foreach define,$(DEFINES),-D$(define))
LD = cc
LFLAGS = -lm -O2
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
         norm.o \
         export.o \
         print.o \
         sparse_matrix.o \
         iterative_method.o \
         functions.o
OFILES = $(foreach fname,$(ONAMES),$(BUILDDIR)/program/$(fname))

%/create-stamp:
	mkdir -p $*
	touch $@

gas: $(OFILES) $(LASPACKOFILES)
	$(LD) -o gas $(OFILES) $(LASPACKOFILES) $(LFLAGS)

$(BUILDDIR)/laspack/%.o: $(BUILDDIR)/laspack/create-stamp lib/laspack/%.c
	$(CC) $(CFLAGS) -o $@ -c lib/laspack/$*.c

$(BUILDDIR)/program/%.o: $(BUILDDIR)/program/create-stamp src/%.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c src/$*.c

clean:
	rm -rf $(BUILDDIR)

distclean: clean
	rm -f gas

.SECONDARY: $(BUILDDIR)/laspack/create-stamp $(BUILDDIR)/program/create-stamp
