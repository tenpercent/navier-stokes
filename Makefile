BUILDDIR = build
DEFINES =
CC = cc
ICC = icc
CFLAGS = -Wall -g -Ilib -Ilib/laspack -std=c99
ICFLAGS = -Wall -g -O3 -scalar-rep -Ilib -Ilib/laspack -std=c99 -mssse3
CPPFLAGS = $(foreach define,$(DEFINES),-D$(define))
LD = cc
ILD = xild
LFLAGS = -lm -g -O2
ILFLAGS = -lm -lirc -limf -lsvml -macosx_version_min 10.9 -L /opt/intel/composer_xe_2013_sp1/lib

# CC = $(ICC)
# CFLAGS = $(ICFLAGS)
# LD = $(ILD)
# LFLAGS = $(ILFLAGS)

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

ifeq (,$(filter NO_LASPACK,$(DEFINES)))
gas: $(OFILES) $(LASPACKOFILES)
	$(LD) -o gas $(OFILES) $(LASPACKOFILES) $(LFLAGS)
else
gas: $(OFILES)
	$(LD) -o gas $(OFILES) $(LFLAGS)
endif

$(BUILDDIR)/laspack/%.o: $(BUILDDIR)/laspack/create-stamp lib/laspack/%.c
	$(CC) $(CFLAGS) -o $@ -c lib/laspack/$*.c

$(BUILDDIR)/program/%.o: $(BUILDDIR)/program/create-stamp src/%.c Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c src/$*.c

clean:
	rm -rf $(BUILDDIR)

distclean: clean
	rm -f gas

.SECONDARY: $(BUILDDIR)/laspack/create-stamp $(BUILDDIR)/program/create-stamp
