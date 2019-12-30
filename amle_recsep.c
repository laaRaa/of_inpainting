/* Copyright (c) 2019 Lara Raad <lara.raad@upf.edu>,
			     <enric.meinhardt@cmla.ens-cachan.fr>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

Smooth Contours is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iio.h"
#include "lib_amle_recsep.h"


// command-line interface
int main(int argc, char *argv[])
{
	if (argc != 11) {
		fprintf(stderr, "usage:\n\t"
		"%s NS lambda err_threshold distance_type neighbourgood_type neighbourhood_ratio data.png mask.png out.png guide.png\n", *argv);
		//0 1  2      3             4             5                  6                   7        8        9       10
                fprintf(stderr, "\n");
                fprintf(stderr, "distance_type = 1, 2, 3 or 4\n");
                fprintf(stderr, "neighbourhood_type = 1 or 2\n");
                fprintf(stderr, "neighbourhood_ratio = 1, 2, 3, 4 or 5\n");
		return 1;
	}
	int niter = 5000; //atoi(argv[1]);
	int nscales = atoi(argv[1]);
	float lambda = atof(argv[2]);
        float err_thresh = atof(argv[3]);
        int dist_type = atoi(argv[4]);
        int nn_type = atoi(argv[5]);
	int r = atoi(argv[6]);
	char *filename_in = argv[7];
	char *filename_mask = argv[8];
	char *filename_out = argv[9];
	char *filename_guide = argv[10];

	// check distance type validity
	if (dist_type<1 || dist_type>4)
		return printf("The distance type should be 1, 2, 3 or 4\n");

	// check neighbourhood type validity
	if (nn_type<1 || nn_type>2)
		return printf("The neighbourhood type should be  1 or 2\n");

	// set neighbourhood size given a valid dx
	// ...
	if (r<1 || r>5)
		return fprintf(stderr, "The neighbourhood ratio should be 1, 2, 3, 4, or 5\n");
        int nn;
        if (nn_type==1)
		switch(r)
		{
			case 1: nn = 8;
				break;
			case 2:
				nn = 16;
				break;
			case 3:
				nn = 32;
				break;
			case 4:
				nn = 56;
				break;
			case 5:
				nn = 88;
				break;
		}
        if (nn_type==2)
		switch(r)
		{
			case 1:
				nn = 8;
				break;
			case 2:
				nn = 24;
				break;
			case 3:
				nn = 48;
				break;
			case 4:
				nn = 80;
				break;
			case 5:
				nn = 120;
				break;
		}

	int w[3], h[3], pd[3];
	float *in = iio_read_image_float_split(filename_in, w, h, pd);
	int w_c, h_c, pd_c;
	float *in_copy = iio_read_image_float_split(filename_in, &w_c, &h_c, &pd_c);

	float *mask = NULL;
	if (0 != strcmp(filename_mask, "NAN")) {
		mask = iio_read_image_float(filename_mask, w+1, h+1);
		if (w[0] != w[1] || h[0] != h[1])
			return fprintf(stderr, "image and mask size mismatch");
	}

	float *guide = iio_read_image_float_split(filename_guide, w+2, h+2, pd+2);
	if (w[0] != w[2] || h[0] != h[2])
		return fprintf(stderr, "image and guide size mismatch");

	float *out = malloc(*w**h**pd*sizeof*out);

	int num_pixels_to_inpaint = 0;
	for (int i = 0; i < *w * *h; i++)
		if (mask && mask[i] > 0)
		{
			num_pixels_to_inpaint += 1;
			for (int l = 0; l < *pd; l++)
				in[*w**h*l+i] = NAN;
		}
	amle_recursive_separable(out, in, guide, *w, *h, pd[0], pd[2], niter, nscales, lambda, err_thresh, dist_type, nn_type, nn);


	// compute EPE between input flow and output flow
	float diff = 0;
	for (int i=0; i<*w**h; i++)
		if (mask[i] > 0)
			diff += sqrt((in_copy[i]-out[i])*(in_copy[i]-out[i]) + (in_copy[i+*w**h]-out[i+*w**h])*(in_copy[i+*w**h]-out[i+*w**h])) ;
	float EPE = diff/num_pixels_to_inpaint;
	printf("%f\n",EPE);

	iio_write_image_float_split(filename_out, out, *w, *h, *pd);

	return 0;
}
