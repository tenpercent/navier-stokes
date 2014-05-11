all: build

build::
	[ -f $@/Makefile ] || (cd $@ && cmake ..)
	$(MAKE) -C $@ $(MAKECMDGOALS) VERBOSE=1

clean: build
	find . -name '*.o' -delete
