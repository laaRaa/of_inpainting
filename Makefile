# compiler and linker options
CFLAGS ?= -O3 -march=native
LDLIBS  = -lpng -ltiff -ljpeg -lm
MINTERP = octave

# variables
OBJ     = amle_recsep.o iio.o
BIN     = amle_recsep

# default target
default: $(BIN)

# build rule
$(BIN): $(OBJ)

# test
test: $(BIN)
	./amle_recsep 3 0.001 0.001 3 1 2 test/input.flo test/mask.png test/out.flo test/guide.png
	$(MINTERP) test_LB.m || true

# bureaucracy
clean:  ; $(RM) $(BIN) $(OBJ) test/out.flo
.PHONY: default test clean


# hack to add -fopenmp only when the CXX compiler supports it
OMPVER := $(shell $(CXX) -dM -E -fopenmp - 2>/dev/null </dev/null |grep _OPENMP)
ifneq ($(OMPVER),)
LDLIBS := $(LDLIBS) -fopenmp
CFLAGS := $(CFLAGS) -fopenmp
endif
