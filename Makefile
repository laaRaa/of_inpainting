CFLAGS ?= -O3 -march=native
LDLIBS  = -lpng -ltiff -ljpeg -lm

BIN = amle_recsep
OBJ = amle_recsep.o iio.o

$(BIN) : $(OBJ)

clean: ; $(RM) $(BIN) $(OBJ)

