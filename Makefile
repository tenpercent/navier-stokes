all: build

build::
	$(MAKE) -C $@ $(MAKECMDGOALS) VERBOSE=1

clean: build
	find . -name '*.o' -delete
