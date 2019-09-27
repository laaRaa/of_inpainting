# compiler and linker options
CFLAGS ?= -O3 -march=native        -fopenmp
LDLIBS  = -lpng -ltiff -ljpeg -lm  -fopenmp

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

# bureaucracy
clean:  ; $(RM) $(BIN) $(OBJ) test/out.flo
.PHONY: default test clean
