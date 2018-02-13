
BIN=.

program = $(BIN)/xreformat_hdf5

include ./config.h

objs = trim_comments_mod.o re_format_hdf5.o

.SUFFIXES: .o .f90

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)

.f90.o:

	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean

clean:
	rm -f $(program) $(objs) *.lst *.o *.mod

all: $(program)
