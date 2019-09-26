CFLAGS ?= -O3 -march=native
LDLIBS  = -lpng -ltiff -ljpeg -lm

<<<<<<< HEAD

SRC    = amle_recsep.c iio.c
BIN    = amle_recsep
COPT   = -std=c99 -O3
OPENMP =

default: $(BIN)

openmp: OPENMP=-fopenmp
openmp: default
=======
BIN = amle_recsep
OBJ = amle_recsep.o iio.o

$(BIN) : $(OBJ)

clean: ; $(RM) $(BIN) $(OBJ)
>>>>>>> f2e7ab5ac8b8e98a0c5007b04164045f5ede326e

$(BIN): $(SRC)
	$(CC) $(COPT) $(OPENMP) -o $@ $^ -lpng -ltiff -ljpeg -lm

clean:
	$(RM) $(BIN)
