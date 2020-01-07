/* 
Modified version: Copyright (c) 2019 Lara Raad <lara.raad@upf.edu>,
Original version: Copyright (c) 2017 Enric Meinhardt-Llopis <enric.meinhardt@cmla.ens-cachan.fr>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

Smooth Contours is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>. 
*/

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "smapa.h"


// utility function that always returns a valid pointer to memory
static void *xmalloc(size_t n)
{
	void *new = malloc(n);
	if (!new)
	{
		fprintf(stderr, "xmalloc: can not malloc %zu bytes\n", n);
		exit(1);
	}
	return new;
}

// the type of a "getpixel" function
typedef float (*getpixel_operator)(float*,int,int,int,int,int);

// extrapolate by 0
inline static float getpixel_0(float *x, int w, int h, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i+j*w+l*w*h];
}

// extrapolate by nearest value (useful for Neumann boundary conditions)
inline static float getpixel_1(float *x, int w, int h, int i, int j, int l)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w+l*w*h];
}


// build a mask of the NAN positions on image "x"
// the output "mask[i][2]" contains the two coordinates of the ith masked pixel
static int (*build_mask(int *out_nmask, float *x, int w, int h))[2]
{
  int nmask = 0;
  for (int i = 0; i < w*h; i++)
    if (isnan(x[i]))
      nmask += 1;
  int (*mask)[2] = xmalloc(w*h*2*sizeof(int)), cx = 0;
  for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++)
      if (isnan(x[j*w + i])) {
        mask[cx][0] = i;
        mask[cx][1] = j;
        cx += 1;
      }
  assert(cx == nmask);

  *out_nmask = nmask;
  return mask;
}

  static int(*get_neighbours(int nn_type, int nn))[2]
{
  int (*n)[2] = xmalloc(nn*2*sizeof(int));

  assert( nn_type==1 || nn_type==2 || nn_type==3 );

  if (nn_type==1)
  {
    int p[][2] = { // {x, y}
      {+1,0}, {0,+1}, {-1,0}, {0,-1},
      {+1,+1}, {-1,-1}, {-1,+1}, {+1,-1}, // 8-neighbors, dx=sqrt(2)
      {+2,+1}, {+1,+2}, {+2,-1}, {+1,-2},
      {-2,-1}, {-1,-2}, {-2,+1}, {-1,+2}, // 16-neighbors, dx = sqrt(3)
      {+3,+1}, {+1,+3}, {+3,-1}, {+1,-3},
      {-3,-1}, {-1,-3}, {-3,+1}, {-1,+3},
      {+3,+2}, {+2,+3}, {+3,-2}, {+2,-3},
      {-3,-2}, {-2,-3}, {-3,+2}, {-2,+3}, // 32-neighbors, dx = sqrt(5)
      {+4,+1}, {+4,+2}, {+4,+3},{+4,-1}, {+4,-2}, {+4,-3},
      {-4,+1}, {-4,+2}, {-4,+3},{-4,-1}, {-4,-2}, {-4,-3},
      {+1,+4}, {+2,+4}, {+3,+4},{-1,+4}, {-2,+4}, {-3,+4},
      {+1,-4}, {+2,-4}, {+3,-4},{-1,-4}, {-2,-4}, {-3,-4}, // 56-neighbors, dx=sqrt(7)
      {+5,+1}, {+5,+2}, {+5,+3},{+5,+4},{+5,-1}, {+5,-2}, {+5,-3},{+5,-4},
      {-5,+1}, {-5,+2}, {-5,+3},{-5,+4},{-5,-1}, {-5,-2}, {-5,-3},{-5,-4},
      {+1,+5}, {+2,+5}, {+3,+5},{+4,+5},{-1,+5}, {-2,+5}, {-3,+5},{-4,+5},
      {+1,-5}, {+2,-5}, {+3,-5},{+4,-5},{-1,-5}, {-2,-5}, {-3,-5},{-4,-5}, // 88-neighbors, dx = 3

    };
    assert( nn >= 0 && nn <= 32 );
    for (int i = 0; i < nn; i++)
    {
      n[i][0] = p[i][0];
      n[i][1] = p[i][1];
    }
  }
  else // nn_type 2 or 3
  {
    int p[][2] = { // {x, y}
      {+1,0}, {0,+1}, {-1,0}, {0,-1},
      {+1,+1}, {-1,-1}, {-1,+1}, {+1,-1}, // 8-neighbors, dx=sqrt(2)
      {+2,+1}, {+1,+2}, {+2,-1}, {+1,-2}, {+2,0}, {+2,+2}, {+2,-2},
      {-2,-1}, {-1,-2}, {-2,+1}, {-1,+2}, {-2,0}, {-2,-2}, {-2,+2},
      {0,-2}, {0,+2}, // 24-neighbors, dx = 2
      {+3,+1}, {+1,+3}, {+3,-1}, {+1,-3},
      {-3,-1}, {-1,-3}, {-3,+1}, {-1,+3},
      {+3,+2}, {+2,+3}, {+3,-2}, {+2,-3},
      {-3,-2}, {-2,-3}, {-3,+2}, {-2,+3},
      {+3,0}, {+3,+3}, {0,+3}, {-3,+3},
      {-3,0}, {-3,-3}, {0,-3}, {+3,-3}, // 48-neighbors, dx = sqrt(6)
      {+4,0}, {+4,+1}, {+4,+2}, {+4,+3}, {+4,+4},
      {-4,0}, {-4,+1}, {-4,+2}, {-4,+3}, {-4,+4},
      {+4,-1}, {+4,-2}, {+4,-3}, {+4,-4},
      {-4,-1}, {-4,-2}, {-4,-3}, {-4,-4},
      {-3,-4}, {-2,-4}, {-1,-4}, {0,-4},
      {+3,-4}, {+2,-4}, {+1,-4},
      {-3,+4}, {-2,+4}, {-1,+4}, {0,+4},
      {+3,+4}, {+2,+4}, {+1,+4}, // 80-neighbors, dx=sqrt(8)
      {+5,0}, {+5,+1}, {+5,+2}, {+5,+3}, {+5,+4}, {+5,+5},
      {-5,0}, {-5,+1}, {-5,+2}, {-5,+3}, {-5,+4}, {-5,+5},
      {+5,-1}, {+5,-2}, {+5,-3}, {+5,-4}, {+5,-5},
      {-5,-1}, {-5,-2}, {-5,-3}, {-5,-4}, {-5,-5},
      {-4,-5}, {-3,-5}, {-2,-5}, {-1,-5}, {0,-5},
      {+4,-5}, {+3,-5}, {+2,-5}, {+1,-5},
      {-4,+5}, {-3,+5}, {-2,+5}, {-1,+5}, {0,-5},
      {+4,+5}, {+3,+5}, {+2,+5}, {+1,+5}, // 120-neighbors, dx=sqrt(10)
    };

    assert( nn >= 0 && nn <= 120 );

    int pi;

    // nn_type = 2 all neighbours are used up to nn
    // nn_type = 3 only nn neighbours are used corresponding to the external square ring
    if (nn_type == 2 || nn==8) pi = 0;
    else if (nn==16) pi = 8;
    else if (nn==24) pi = 24;
    else if (nn==32) pi = 48;
    else if (nn==40) pi = 80;

    for (int i = pi; i < pi+nn; i++)
    {
      n[i-pi][0] = p[i][0];
      n[i-pi][1] = p[i][1];
    }
  }

  return n;
}

// evaluate an image in a neighborhood
static int get_nvals(float *v, float *wv2, float *x, float *weight, int w, int h, int i, int j,
    int nn_type, int nn)
{

  int r = 0;
  int (*n)[2] = get_neighbours(nn_type, nn);
  int (*pn)[2] = n;

  for (int p = 0; p < nn; p++)
  {
    int ii = i + pn[p][0];
    int jj = j + pn[p][1];
    if (ii >= 0 && jj >= 0 && ii < w && jj < h)
    {
      v[r] = x[w*jj+ii];
      if (wv2)
        wv2[r] = weight[j*w+i+h*w*p];
      r += 1;
    }
  }

  free(n);

  return r;
}

static void get_eikonal_idx(int *ind_pos, int *ind_neg, float *x, float *w, int n, float x0)
{
	*ind_pos = 0;
        *ind_neg = 0;

        float eik_pos_max = -INFINITY;
        float eik_neg_min = +INFINITY;

        for (int j=0; j<n; j++)
        {
          float eik = ((x[j]-x0))*(w[j]);
          if (eik > eik_pos_max)
          {
            eik_pos_max = eik;
            *ind_pos = j;
          }
          if (eik < eik_neg_min)
          {
            eik_neg_min = eik;
            *ind_neg = j;
          }
        }
}

static float amle_iteration(float *x, float *weight, int w, int h, int (*mask)[2], int nmask, int nn_type, int nn)
{
	float actus = 0;
	float actumax = 0;

#ifdef _OPENMP
#pragma omp parallel for
#endif//_OPENMP
	for (int p = 0; p < nmask; p++)
	{
		int i = mask[p][0];
		int j = mask[p][1];
		int idx = j*w + i, indy, indz;
                float x0 = x[idx];
		float value[0x100], weights[0x100];
		int nv = get_nvals(value, weights, x, weight, w, h, i, j, nn_type, nn);
                get_eikonal_idx(&indy, &indz, value, weights, nv, x0);
		float a = weights[indz];
		float b = weights[indy];
		float newx = (a*value[indz] + b*value[indy]) / (a + b);

		if (fabs(x[idx]-newx) > actumax)
			actumax = fabs(x[idx]-newx);
		x[idx] = newx;
	}
	return actumax;
}

// fill the holes of the image x using an infinity harmonic function
static void amle_extension(
		float *y,        // output image
		float *x,        // input image (NAN values indicate holes)
		float *weight,
		int w,           // image width
		int h,           // image height
		int niter,       // number of iterations to run
		float *initialization,
                float err_thresh,
                int nn_type,
		int nn,
                int scale,       // current scale
                int scale_num    // total number of scales
		)
{
	// build list of masked pixels
	int nmask, (*mask)[2] = build_mask(&nmask, x, w, h);

	// initialize result of u at it (k-1)
	float *y_old = xmalloc(w*h*sizeof*y_old);
	
	// initialize the solution to the given data at the masked pixels
	for (int i = 0; i < w*h; i++)
	{
		y[i] = isfinite(x[i]) ? x[i] : initialization[i];
		y_old[i] = isfinite(x[i]) ? x[i] : initialization[i];
	}

        int count = 0;
        float e_k = 1000;

	// do the requested iterations
	while (count < niter && e_k > err_thresh)
	{
		float u = amle_iteration(y, weight, w, h, mask, nmask, nn_type, nn);
		float diff = 0;
		for (int i=0; i<w*h; i++)
			if (!isfinite(x[i]))
				diff += fabs(y[i]-y_old[i]);
		e_k = diff/nmask;
		for (int i=0; i<w*h; i++)
			y_old[i] = y[i];
		count += 1;
	}

	printf("niter %d, e_k %g\n", count, e_k);

	free(mask);
	free(y_old);
}

// zoom-out by 2x2 block averages
// NANs are discarded when possible
static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih, int pd)
{
	getpixel_operator p = getpixel_1;
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
        for (int l=0; l<pd; l++)
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	{
		float a[4], m = 0;
		a[0] = p(in, iw, ih, 2*i, 2*j, l);
		a[1] = p(in, iw, ih, 2*i+1, 2*j, l);
		a[2] = p(in, iw, ih, 2*i, 2*j+1, l);
		a[3] = p(in, iw, ih, 2*i+1, 2*j+1, l);
		int cx = 0;
		for (int k = 0; k < 4; k++)
			if (isfinite(a[k])) {
				m += a[k];
				cx += 1;
			}
		out[ow*j + i + l*ow*oh] = cx ? m/cx : NAN;
	}
}

// bilinear_interpolation
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}
static float bilinear_interpolation(float *x, int w, int h, float p, float q, int l)
{
	int ip = p;
	int iq = q;
	float a = getpixel_1(x, w, h, ip, iq, l);
	float b = getpixel_1(x, w, h, ip+1, iq, l);
	float c = getpixel_1(x, w, h, ip  , iq+1, l);
	float d = getpixel_1(x, w, h, ip+1, iq+1, l);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

// zoom-in by replicating pixels into 2x2 blocks
// no NAN's are expected in the input image
static void zoom_in_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih, int pd)
{
	//getpixel_operator p = getpixel_1;
	assert(abs(2*iw-ow) < 2);
	assert(abs(2*ih-oh) < 2);
        for (int l = 0; l < pd; l++)
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
		//out[ow*j+i + l*ow*oh] = p(in, iw, ih, round((i-0.5)/2), round((j-0.5)/2), l);
                out[ow*j+i] = bilinear_interpolation(in, iw, ih, (i-0.5)/2.0, (j-0.5)/2.0, l);
}

// compute the distances of each point of the image to all
// neighbours. The distance is
// d(x,y) = sqrt((1-lambda)*|I(x)-I(y)|^2 + lambda*|x-y|^2)
static void compute_weight_1(float *weight, float *u, int w, int h, int pd, float lambda, int nn_type, int nn)
{

  int (*n)[2] = get_neighbours(nn_type,nn);

  for(int j=0; j<h; j++)
    for(int i=0; i<w; i++)
      for (int p = 0; p < nn; p++)
      {
        int ii = i + n[p][0];
        int jj = j + n[p][1];
        if (ii >= 0 && jj >= 0 && ii < w && jj < h)
        {
          weight[j*w+i+h*w*p] = 0;
	  float a = 0;
	  float b = (ii-i)*(ii-i) + (jj-j)*(jj-j);
          for (int l=0; l<pd; l++)
            a += (u[jj*w+ii+w*h*l]-u[j*w+i+w*h*l])*(u[jj*w+ii+w*h*l]-u[j*w+i+w*h*l]);
	  a /= pd;
          weight[j*w+i+h*w*p] = 1/(sqrt((1-lambda)*a + lambda*b));
        }
        if (ii < 0 || jj < 0 || ii >= w || jj >= h)
          weight[j*w+i+h*w*p] = -1;
      }

  free(n);
}

// compute the distances of each point of the image to all
// neighbours. The distance is
// d(x,y) = (1-lambda)*|I(x)-I(y)| + lambda*|x-y|
static void compute_weight_2(float *weight, float *u, int w, int h, int pd, float lambda, int nn_type, int nn)
{

  int (*n)[2] = get_neighbours(nn_type,nn);

  for(int j=0; j<h; j++)
    for(int i=0; i<w; i++)
      for (int p = 0; p < nn; p++)
      {
        int ii = i + n[p][0];
        int jj = j + n[p][1];
        if (ii >= 0 && jj >= 0 && ii < w && jj < h)
        {
	  float a = 0;
	  float b = sqrt((ii-i)*(ii-i) + (jj-j)*(jj-j));
          for (int l=0; l<pd; l++)
            a += (u[jj*w+ii+w*h*l]-u[j*w+i+w*h*l])*(u[jj*w+ii+w*h*l]-u[j*w+i+w*h*l]);
	  a /= pd;
	  a = sqrt(a);
          weight[j*w+i+h*w*p] = 1/((1-lambda)*a + lambda*b);
        }
        if (ii < 0 || jj < 0 || ii >= w || jj >= h)
          weight[j*w+i+h*w*p] = -1;
      }

  free(n);
}

// compute the distances of each point of the image to all
// neighbours. The distance is
// d(x,y) = (1-lambda)*|I(x)-I(y)|^2 + lambda*|x-y|^2
static void compute_weight_3(float *weight, float *u, int w, int h, int pd, float lambda, int nn_type, int nn)
{

  int (*n)[2] = get_neighbours(nn_type,nn);

  for(int j=0; j<h; j++)
    for(int i=0; i<w; i++)
      for (int p = 0; p < nn; p++)
      {
        int ii = i + n[p][0];
        int jj = j + n[p][1];
        if (ii >= 0 && jj >= 0 && ii < w && jj < h)
        {
	  float a = 0;
	  float b = (ii-i)*(ii-i) + (jj-j)*(jj-j);
          for (int l=0; l<pd; l++)
            a += (u[jj*w+ii+w*h*l]-u[j*w+i+w*h*l])*(u[jj*w+ii+w*h*l]-u[j*w+i+w*h*l]);
	  a /= pd;
	  weight[j*w+i+h*w*p] = 1/((1-lambda)*a + lambda*b);
        }
        if (ii < 0 || jj < 0 || ii >= w || jj >= h)
          weight[j*w+i+h*w*p] = -1;
      }

  free(n);
}

// compute the distances of each point of the image to all
// neighbours. The distance is
// d(x,y) = (1-lambda)*||patch(x)-patch(y)||^2 + lambda*||x-y||^2
// where patch(x) is the patch centered at x of size sxs
// s has to be odd
static void compute_weight_4(float *weight, float *u, int w, int h, int pd, float lambda, int nn_type, int nn)
{

  int s = 5; // add as parameter to the function: int s
  int (*n)[2] = get_neighbours(nn_type,nn);

  for(int j=0; j<h; j++)
    for(int i=0; i<w; i++)
    {
      for (int p = 0; p < nn; p++)
      {
        int ii = i + n[p][0];
        int jj = j + n[p][1];
	if (ii >= 0 && jj >= 0 && ii < w && jj < h)
        {
          float a = 0;
          for (int l=0; l<pd; l++)
            for (int m=-(s-1)/2; m<= (s-1)/2; m++)
	      for (int n=-(s-1)/2; n<= (s-1)/2; n++)
		if (jj+n>=0 && jj+n<h && ii+m>=0 && ii+m<w && j+n>=0 && j+n<h && i+m>=0 && i+m<w)
		  a += (u[(jj+n)*w+(ii+m)+w*h*l]-u[(j+n)*w+(i+m)+w*h*l])*(u[(jj+n)*w+(ii+m)+w*h*l]-u[(j+n)*w+(i+m)+w*h*l]);
          a = a/s/s/pd; // normalize by the number of pixels in a patch
          float b = (ii-i)*(ii-i) + (jj-j)*(jj-j);
          weight[j*w+i+h*w*p] = 1/((1-lambda)*a + lambda*b);
        }
        if (ii < 0 || jj < 0 || ii >= w || jj >= h)
          weight[j*w+i+h*w*p] = -1;
      }
    }

  free(n);
}


static void compute_weight_map(float *weight, float *u, int w, int h, int pd, float lambda, int w_type, int nn_type, int nn)
{
  switch (w_type)
  {
	case 1:
	  compute_weight_1(weight, u, w, h, pd, lambda, nn_type, nn);
	  break;
	case 2:
	  compute_weight_2(weight, u, w, h, pd, lambda, nn_type, nn);
	  break;
	case 3:
	  compute_weight_3(weight, u, w, h, pd, lambda, nn_type, nn);
	  break;
	case 4:
	  compute_weight_4(weight, u, w, h, pd, lambda, nn_type, nn);
	  break;
  }
}
#include "iio.h"

SMART_PARAMETER_SILENT(AMLE_ONLY,0)

void amle_recursive(float *out, float *in, float *guide, int w, int h, int niter, int scale, float lambda, int pd_guide, float err_thresh, int w_type, int nn_type, int nn, int scale_num)
{
  float *weight = xmalloc(w*h*nn*sizeof*weight);
  compute_weight_map(weight, guide, w, h, pd_guide, lambda, w_type, nn_type, nn);
  float *init = xmalloc(w*h*sizeof*init);
  if (scale > 1)
  {
    int ws = ceil(w/2.0);
    int hs = ceil(h/2.0);
    float *ins = xmalloc(ws*hs*sizeof*ins);
    float *guides = xmalloc(ws*hs*pd_guide*sizeof*guides);
    float *outs = xmalloc(ws*hs*sizeof*outs);
    zoom_out_by_factor_two(ins, ws, hs, in, w, h, 1);
    zoom_out_by_factor_two(guides, ws, hs, guide, w, h, pd_guide);
    amle_recursive(outs, ins, guides, ws, hs, niter, scale - 1, lambda, pd_guide, err_thresh, w_type, nn_type, nn, scale_num);
    zoom_in_by_factor_two(init, w, h, outs, ws, hs, 1);

    free(ins);
    free(outs);
    free(guides);
  } else {
    for (int i = 0 ; i < w*h; i++)
      init[i] = 0;
  }

  if (AMLE_ONLY() > 0 && AMLE_ONLY()!=w) niter = 0;
  amle_extension(out, in, weight, w, h, niter, init, err_thresh, nn_type, nn, scale, scale_num);
  free(init);
  free(weight);
}


// extension by AMLE each channel of a color image
void amle_recursive_separable(float *out, float *in, float *guide, int w, int h, int pd, int pd_guide, int niter, int nscale, float lambda, float err_thresh, int w_type, int nn_type, int nn)
{
	for (int l = 0; l < pd; l++)
	{
		float *outl = out + w*h*l;
		float *inl = in + w*h*l;
		amle_recursive(outl, inl, guide, w, h, niter, nscale, lambda, pd_guide, err_thresh, w_type, nn_type, nn, nscale);
	}
}
