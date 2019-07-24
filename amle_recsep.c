#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN(a,b) ( (a)<(b) ? (a) : (b) )
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )
#define ABS(a) (a<0)?-(a):(a)

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

#include "smapa.h"
SMART_PARAMETER_SILENT(AMLE_NN,8)

  static int(*get_neighbours(int nn_type))[2]
{
  int nn = AMLE_NN();
  int (*n)[2] = xmalloc(nn*2*sizeof(int));

  assert( nn_type==1 || nn_type==2 || nn_type==3 );

  if (nn_type==1)
  {
    int p[][2] = { // {x, y}
      {+1,0}, {0,+1}, {-1,0}, {0,-1}, // 4-connexity
      {+1,+1}, {-1,-1}, // 6-connexity
      {-1,+1}, {+1,-1}, // 8-connexity
      {+2,+1}, {+1,+2}, {+2,-1}, {+1,-2},
      {-2,-1}, {-1,-2}, {-2,+1}, {-1,+2}, // 16-neighbors
      {+3,+1}, {+1,+3}, {+3,-1}, {+1,-3},
      {-3,-1}, {-1,-3}, {-3,+1}, {-1,+3}, // 24-neighbors
      {+3,+2}, {+2,+3}, {+3,-2}, {+2,-3},
      {-3,-2}, {-2,-3}, {-3,+2}, {-2,+3}, // 32-neighbors
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
      {+1,0}, {0,+1}, {-1,0}, {0,-1}, // 4-connexity
      {+1,+1}, {-1,-1}, // 6-connexity
      {-1,+1}, {+1,-1}, // 8-connexity
      {+2,+1}, {+1,+2}, {+2,-1}, {+1,-2}, {+2,0}, {+2,+2}, {+2,-2},
      {-2,-1}, {-1,-2}, {-2,+1}, {-1,+2}, {-2,0}, {-2,-2}, {-2,+2},
      {0,-2}, {0,+2}, // 24-neighbours
      {+3,+1}, {+1,+3}, {+3,-1}, {+1,-3},
      {-3,-1}, {-1,-3}, {-3,+1}, {-1,+3},
      {+3,+2}, {+2,+3}, {+3,-2}, {+2,-3},
      {-3,-2}, {-2,-3}, {-3,+2}, {-2,+3},
      {+3,0}, {+3,+3}, {0,+3}, {-3,+3},
      {-3,0}, {-3,-3}, {0,-3}, {+3,-3}, // 48-neighbors
      {+4,0}, {+4,+1}, {+4,+2}, {+4,+3}, {+4,+4},
      {-4,0}, {-4,+1}, {-4,+2}, {-4,+3}, {-4,+4},
      {+4,-1}, {+4,-2}, {+4,-3}, {+4,-4},
      {-4,-1}, {-4,-2}, {-4,-3}, {-4,-4},
      {-3,-4}, {-2,-4}, {-1,-4}, {0,-4},
      {+3,-4}, {+2,-4}, {+1,-4},
      {-3,+4}, {-2,+4}, {-1,+4}, {0,+4},
      {+3,+4}, {+2,+4}, {+1,+4}, // 80-neighbours
      {+5,0}, {+5,+1}, {+5,+2}, {+5,+3}, {+5,+4}, {+5,+5},
      {-5,0}, {-5,+1}, {-5,+2}, {-5,+3}, {-5,+4}, {-5,+5},
      {+5,-1}, {+5,-2}, {+5,-3}, {+5,-4}, {+5,-5},
      {-5,-1}, {-5,-2}, {-5,-3}, {-5,-4}, {-5,-5},
      {-4,-5}, {-3,-5}, {-2,-5}, {-1,-5}, {0,-5},
      {+4,-5}, {+3,-5}, {+2,-5}, {+1,-5},
      {-4,+5}, {-3,+5}, {-2,+5}, {-1,+5}, {0,-5},
      {+4,+5}, {+3,+5}, {+2,+5}, {+1,+5}, // 120-neighbours
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
static int get_nvals(float *v, float *wv2, float *x, float *dist, int w, int h, int i, int j,
    int nn_type)
{

  int r = 0;
  int (*n)[2] = get_neighbours(nn_type);
  int (*pn)[2] = n;
  int nn = AMLE_NN();

  for (int p = 0; p < nn; p++)
  {
    int ii = i + pn[p][0];
    int jj = j + pn[p][1];
    if (ii >= 0 && jj >= 0 && ii < w && jj < h)
    {
      v[r] = x[w*jj+ii];
      if (wv2)
        wv2[r] = dist[j*w+i+h*w*p];
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
          float eik = ((x[j]-x0))/(w[j]);
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

static float amle_iteration(float *x, float *dist, int w, int h, int (*mask)[2], int nmask, float *error, int nn_type)
{
	float actus = 0;
	float actumax = 0;
        float err_aux = 0.0;
        float norm_old = 0.0;
	for (int i = 0; i < w*h; i++)
          norm_old += x[i]*x[i];

        error[0] = 0; // initialize error to 0
//#pragma omp parallel for
	for (int p = 0; p < nmask; p++)
	{
		int i = mask[p][0];
		int j = mask[p][1];
		int idx = j*w + i, indy, indz;
                float x0 = x[idx];
		float value[0x100], weight[0x100];
		int nv = get_nvals(value, weight, x, dist, w, h, i, j, nn_type);
                get_eikonal_idx(&indy, &indz, value, weight, nv, x0);
		float a = weight[indz];
		float b = weight[indy];
		float newx = (a*value[indy] + b*value[indz]) / (a + b);

                err_aux = fabs(newx-x[idx]);
                if (err_aux>error[0])
                  error[0] = err_aux;

		actus += fabs(x[idx] - newx);
		if (fabs(x[idx]-newx) > actumax)
			actumax = fabs(x[idx]-newx);
		x[idx] = newx;
	}
	return actumax;
}

// fill the holes of the image x using an infinity harmonic function
static void inf_harmonic_extension_with_init(
		float *y,        // output image
		float *x,        // input image (NAN values indicate holes)
		float *dist,
		int w,           // image width
		int h,           // image height
		int niter,       // number of iterations to run
		float *initialization,
                float err_thresh,
                int nn_type,
                int scale,       // current scale
                int scale_num    // total number of scales
		)
{
	// build list of masked pixels
	int nmask, (*mask)[2] = build_mask(&nmask, x, w, h);
        float err_thresh_aux = err_thresh; // relax threshold value for convergence in coarser scales

	// initialize the solution to the given data at the masked pixels
	for (int i = 0; i < w*h; i++)
          y[i] = isfinite(x[i]) ? x[i] : initialization[i];

        int counter = 0;
        float err[] = {1000};

        // do the requested iterations
        if (scale==scale_num)
          /* while (counter < niter && err[0] > err_thresh) */
          while (counter < niter)
          {
            float u = amle_iteration(y, dist, w, h, mask, nmask, err, nn_type);
            counter += 1;
          }
        else
          while (counter < niter)
          {
            float u = amle_iteration(y, dist, w, h, mask, nmask, err, nn_type);
            counter += 1;
          }

        /* printf("number of iterations %d\n", counter); */
        /* printf("min error between 2 iterations: %g\n", err[0]); */

	free(mask);
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
	getpixel_operator p = getpixel_1;
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
// d(x,y) = (1-lambda)*|I(x)-I(y)|^2 + lambda*|x-y|^2
static void distance0(float *dist, float *u, int w, int h, int pd, float lambda, int nn_type)
{

  int nn = AMLE_NN();
  int (*n)[2] = get_neighbours(nn_type);

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
	  dist[j*w+i+h*w*p] = (1-lambda)*a + lambda*b;
        }
        if (ii < 0 || jj < 0 || ii >= w || jj >= h)
          dist[j*w+i+h*w*p] = -1;
      }

  free(n);
}

// compute the distances of each point of the image to all
// neighbours. The distance is
// d(x,y) = sqrt((1-lambda)*|I(x)-I(y)|^2 + lambda*|x-y|^2)
static void distance1(float *dist, float *u, int w, int h, int pd, float lambda, int nn_type)
{

  int nn = AMLE_NN();
  int (*n)[2] = get_neighbours(nn_type);

  for(int j=0; j<h; j++)
    for(int i=0; i<w; i++)
      for (int p = 0; p < nn; p++)
      {
        int ii = i + n[p][0];
        int jj = j + n[p][1];
        if (ii >= 0 && jj >= 0 && ii < w && jj < h)
        {
          dist[j*w+i+h*w*p] = 0;
	  float a = 0;
	  float b = (ii-i)*(ii-i) + (jj-j)*(jj-j);
          for (int l=0; l<pd; l++)
            a += (u[jj*w+ii+w*h*l]-u[j*w+i+w*h*l])*(u[jj*w+ii+w*h*l]-u[j*w+i+w*h*l]);
	  a /= pd;
          dist[j*w+i+h*w*p] = sqrt((1-lambda)*a + lambda*b);
        }
        if (ii < 0 || jj < 0 || ii >= w || jj >= h)
          dist[j*w+i+h*w*p] = -1;
      }

  free(n);
}

// compute the distances of each point of the image to all
// neighbours. The distance is
// d(x,y) = (1-lambda)*|I(x)-I(y)| + lambda*|x-y|
static void distance2(float *dist, float *u, int w, int h, int pd, float lambda, int nn_type)
{

  int nn = AMLE_NN();
  int (*n)[2] = get_neighbours(nn_type);

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
          dist[j*w+i+h*w*p] = (1-lambda)*a + lambda*b;
        }
        if (ii < 0 || jj < 0 || ii >= w || jj >= h)
          dist[j*w+i+h*w*p] = -1;
      }

  free(n);
}

// compute the distances of each point of the image to all
// neighbours. The distance is
// d(x,y) = (1-lambda)*||patch(x)-patch(y)||^2 + lambda*||x-y||^2
// where patch(x) is the patch centered at x of size sxs
// s has to be odd
static void distance3(float *dist, float *u, int w, int h, int pd, float lambda, int nn_type)
{

  int s = 5; // add as parameter to the function: int s
  int nn = AMLE_NN();
  int (*n)[2] = get_neighbours(nn_type);

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
          a = a/s/s/3; // normalize by the number of pixels in a patch
          float b = (ii-i)*(ii-i) + (jj-j)*(jj-j);
          dist[j*w+i+h*w*p] = (1-lambda)*a + lambda*b;
        }
        if (ii < 0 || jj < 0 || ii >= w || jj >= h)
          dist[j*w+i+h*w*p] = -1;
      }
    }

  free(n);
}


static void distance(float *dist, float *u, int w, int h, int pd, float lambda, int dist_type, int nn_type)
{
  if (dist_type==0)
    distance0(dist, u, w, h, pd, lambda, nn_type);
  else if (dist_type==1)
    distance1(dist, u, w, h, pd, lambda, nn_type);
  else if (dist_type==2)
    distance2(dist, u, w, h, pd, lambda, nn_type);
  else if (dist_type==3)
    distance3(dist, u, w, h, pd, lambda, nn_type);

}
#include "iio.h"

SMART_PARAMETER_SILENT(AMLE_ONLY,0)

void amle_recursive(float *out, float *in, float *guide, float *dist, int w, int h, int niter, int scale, float lambda, int pd_guide, float err_thresh, int dist_type, int nn_type, int scale_num)
{
  float *init = xmalloc(w*h*sizeof*init);
  if (scale > 1)
  {
    int nn = AMLE_NN();
    int ws = ceil(w/2.0);
    int hs = ceil(h/2.0);
    float *ins = xmalloc(ws*hs*sizeof*ins);
    float *guides = xmalloc(ws*hs*pd_guide*sizeof*guides);
    float *dists = xmalloc(ws*hs*nn*sizeof*dists);
    float *outs = xmalloc(ws*hs*sizeof*outs);
    zoom_out_by_factor_two(ins, ws, hs, in, w, h, 1);
    zoom_out_by_factor_two(guides, ws, hs, guide, w, h, pd_guide);
    distance(dists, guides, ws, hs, pd_guide, lambda, dist_type, nn_type);
    amle_recursive(outs, ins, guides, dists, ws, hs, niter, scale - 1, lambda, pd_guide, err_thresh, dist_type, nn_type, scale_num);
    zoom_in_by_factor_two(init, w, h, outs, ws, hs, 1);

    free(ins);
    free(outs);
    free(guides);
    free(dists);
  } else {
    for (int i = 0 ; i < w*h; i++)
      init[i] = 0;
  }

  if (AMLE_ONLY() > 0 && AMLE_ONLY()!=w) niter = 0;
  inf_harmonic_extension_with_init(out, in, dist, w, h, niter, init, err_thresh, nn_type, scale, scale_num);
  free(init);
}


// extension by AMLE each channel of a color image
void amle_recursive_separable(float *out, float *in, float *guide, int w, int h, int pd, int pd_guide, int niter, int nscale, float lambda,
    float err_thresh, int dist_type, int nn_type)
{

	int nn = AMLE_NN();
	float *dist = xmalloc(w*h*nn*sizeof*dist);

	distance(dist, guide, w, h, pd_guide, lambda, dist_type, nn_type);// do this in amle_recursive!!!

	for (int l = 0; l < pd; l++)
	{
		float *outl = out + w*h*l;
		float *inl = in + w*h*l;
		amle_recursive(outl, inl, guide, dist, w, h, niter, nscale, lambda, pd_guide, err_thresh, dist_type, nn_type, nscale);
	}

	free(dist);
}

int main(int argc, char *argv[])
{
	if (argc != 11) {
		fprintf(stderr, "usage:\n\t"
		"%s NITER NS lambda err_threshold distance_type neighbourgood_type data.png mask.png out.png guide.png\n", *argv);
		//0 1     2  3      4             5             6                  7        8        9       10
                fprintf(stderr, "\n");
                fprintf(stderr, "neighbourhood_type = 1 --> AMLE_NN = 4,8,16,24,32\n");
                fprintf(stderr, "neighbourhood_type = 2 --> AMLE_NN = 4,8,24,48,80,120\n");
                fprintf(stderr, "neighbourhood_type = 3 --> AMLE_NN = 8,16,24,40,48\n");
		return 1;
	}
	int niter = atoi(argv[1]);
	int nscales = atoi(argv[2]);
	float lambda = atof(argv[3]);
        float err_thresh = atof(argv[4]);
        int dist_type = atoi(argv[5]);
        int nn_type = atoi(argv[6]);
	char *filename_in = argv[7];
	char *filename_mask = argv[8];
	char *filename_out = argv[9];
	char *filename_guide = argv[10];

        // check distance type validity
        if (dist_type<0 | dist_type>3)
          return printf("The distance type shoulb be 0, 1, 2 or 3\n");

        // check neighbourhood type validity
        if (nn_type<1 | nn_type>3)
          return printf("The neighbourhood type shoulb be  1, 2 or 3\n");

        int nn = AMLE_NN();
        if (nn_type==1)
          if (nn!=4 & nn!=8 & nn!=16 & nn!=24 & nn!=32)
            return fprintf(stderr, "neighbourhood_type = 1 --> AMLE_NN = 4,8,16,24,32\n");
        if (nn_type==2)
          if (nn!=4 & nn!=8 & nn!=24 & nn!=48 & nn!=80 & nn!=120)
            return fprintf(stderr, "neighbourhood_type = 2 --> AMLE_NN = 4,8,24,48,80,120\n");
        if (nn_type==3)
          if (nn!=8 & nn!=16 & nn!=24 & nn!=40 & nn!=48)
            return fprintf(stderr, "neighbourhood_type = 3 --> AMLE_NN = 8,16,24,40,48\n");

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

	float *out = xmalloc(*w**h**pd*sizeof*out);

	int num_pixels_to_inpaint = 0;
	for (int i = 0; i < *w * *h; i++)
		if (mask && mask[i] > 0)
		{
			num_pixels_to_inpaint += 1;
			for (int l = 0; l < *pd; l++)
				in[*w**h*l+i] = NAN;
		}
	amle_recursive_separable(out, in, guide, *w, *h, *pd, pd[2], niter, nscales, lambda, err_thresh, dist_type, nn_type);

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
