/* 
Original version: Copyright (c) 2019 Enric Meinhardt-Llopis <enric.meinhardt@cmla.ens-cachan.fr>

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
