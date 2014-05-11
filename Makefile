CMAKEFLAGS = -DCMAKE_BUILD_TYPE=DEBUG -DPDESolver_COLORFUL_OUTPUT=1

all: build

build::
	[ -f $@/Makefile ] || (cd $@ && cmake .. $(CMAKEFLAGS))
	$(MAKE) -C $@ $(MAKECMDGOALS) VERBOSE=1

clean: build
