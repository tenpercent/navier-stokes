BUILDDIR = build
DEFINES =
CC = cc
CFLAGS = -Wall -g -O2
CPPFLAGS_LASPACK = -Ilib 
CPPFLAGS_PROGRAM = -Ilib $(foreach define,$(DEFINES),-D$(define))
LD = cc
LFLAGS = -lm -g -O2

COMMONONAMES   = calculation.o construction.o initialize.o norm.o export.o \
                 print.o sparse_matrix.o iterative_method.o functions.o
COMMONOFILES   = $(foreach fname,$(COMMONONAMES),$(BUILDDIR)/program/$(fname))
SMOOTHONAMES   = fill_system_smooth.o pdesolver_smooth.o start_conditions_12.o
SMOOTHOFILES   = $(foreach fname,$(SMOOTHONAMES),$(BUILDDIR)/program/$(fname))
NSMOOTH3ONAMES = fill_system_non_smooth.o pdesolver_non_smooth.o start_conditions_3.o
NSMOOTH3OFILES = $(foreach fname,$(NSMOOTH3ONAMES),$(BUILDDIR)/program/$(fname))
NSMOOTH4ONAMES = fill_system_non_smooth.o pdesolver_non_smooth.o start_conditions_4.o
NSMOOTH4OFILES = $(foreach fname,$(NSMOOTH4ONAMES),$(BUILDDIR)/program/$(fname))

ifeq (,$(filter NO_LASPACK,$(DEFINES)))
LASPACKONAMES  = eigenval.o errhandl.o factor.o itersolv.o matrix.o \
                 operats.o precond.o qmatrix.o rtc.o vector.o
LASPACKOFILES  = $(foreach fname,$(LASPACKONAMES),$(BUILDDIR)/laspack/$(fname))
else
LASPACKOFILES  =
endif

all: gas_smooth gas_nonsmooth_3 gas_nonsmooth_4

%/create-stamp:
	mkdir -p $*
	touch $@

gas_smooth: $(COMMONOFILES) $(SMOOTHOFILES) $(LASPACKOFILES)
	$(LD) -o $@ $(COMMONOFILES) $(SMOOTHOFILES) $(LASPACKOFILES) $(LFLAGS)

gas_nonsmooth_3: $(COMMONOFILES) $(NSMOOTH3OFILES) $(LASPACKOFILES)
	$(LD) -o $@ $(COMMONOFILES) $(NSMOOTH3OFILES) $(LASPACKOFILES) $(LFLAGS)

gas_nonsmooth_4: $(COMMONOFILES) $(NSMOOTH4OFILES) $(LASPACKOFILES)
	$(LD) -o $@ $(COMMONOFILES) $(NSMOOTH4OFILES) $(LASPACKOFILES) $(LFLAGS)

$(BUILDDIR)/laspack/%.o: lib/laspack/%.c $(BUILDDIR)/laspack/create-stamp Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS_LASPACK) -o $@ -c $<

$(BUILDDIR)/program/%.o: src/%.c $(BUILDDIR)/program/create-stamp Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS_PROGRAM) -o $@ -c $<

$(BUILDDIR)/tests/%.o: tests/%.c $(BUILDDIR)/tests/create-stamp Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS_PROGRAM) -Isrc -o $@ -c $<

test: $(BUILDDIR)/tests/test_sparse_matrix.o $(BUILDDIR)/program/sparse_matrix.o $(LASPACKOFILES)
	$(LD) -o $(BUILDDIR)/test_sparse_matrix $< $(BUILDDIR)/program/sparse_matrix.o $(LASPACKOFILES) $(LFLAGS)
	./$(BUILDDIR)/test_sparse_matrix && echo "Test run successfully."

pngs:
	./generate_pngs.py results/dat

clean:
	rm -rf $(BUILDDIR)

distclean: clean
	rm -rf results
	rm -f gas_smooth gas_nonsmooth_3 gas_nonsmooth_4

.SECONDARY: $(BUILDDIR)/laspack/create-stamp $(BUILDDIR)/program/create-stamp $(BUILDDIR)/tests/create-stamp
