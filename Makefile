
amle_recsep: amle_recsep.c iio.c iio.h smapa.h
	cc -O3 -o amle_recsep amle_recsep.c iio.c -lpng -ltiff -ljpeg -lm

clean:
	rm -f amle_recsep output.flo

