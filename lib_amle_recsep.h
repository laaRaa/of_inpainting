

// extension by AMLE each channel of a color image
void amle_recursive_separable(
		float *out,        // pre-allocated pointer to output data
		float *in,         // input image with NANs to be filled-in
		float *guide,      // input guide image
		int w,             // width of images
		int h,             // height of images
		int pd,            // number of color channels
		int pd_guide,
		int niter,         // number of AMLE iterations at each scale
		int nscale,        // number of pyramid levels
		float lambda,      // lambda parameter
		float err_thresh,  // error threshold to stop itertions
		int dist_type,     // distance type
		int nn_type,       // neighborhood type
		int nn             // number of neighbors
		);
