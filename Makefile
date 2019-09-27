CFLAGS ?= -O3 -march=native
LDLIBS  = -lpng -ltiff -ljpeg -lm

SRC    = amle_recsep.c iio.c
BIN    = amle_recsep
COPT   = -std=c99 -O3
OPENMP =

default: $(BIN)

openmp: OPENMP=-fopenmp
openmp: default

$(BIN): $(SRC)
	$(CC) $(COPT) $(OPENMP) -o $@ $^ -lpng -ltiff -ljpeg -lm

clean:
	$(RM) $(BIN)
