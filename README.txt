On anisotropic optical flow inpainting algorithms
Version 1.0 26/09/2019
Lara Raad (lara.raad@upf.edu)

README file for two optical flow inpainting algorithms using the AMLE and
Laplace-Beltrami operators on a 2D Riemannian manifold.

-------------------------------
AMLE OPTICAL FLOW INTERPOLATION
-------------------------------

 - Compilation:

   make           : normal compilation
   make openmp    : compilation with OpenMP parallelism support

 - Required libraries:

libpng, libtiff, libjpeg

 - Test:

./amle_recsep 3 0.001 0.001 3 1 2 test/input.flo test/mask.png test/output_amle.flo test/guide.png

 - Usage:

./amle_recsep S lambda epsilon g nt r input.flo mask.png output.flo guide.png

Required parameters:
   S          :   number of scales
   lambda     :   anisotropic weight
   epsilon    :   stopping criterion threshold
   g          :   weight selection
   nt         : type of local neighbourhood
   r          : neighbourhood ratio
   input.flo  : input flow to inpaint
   mask.png   : inpainting mask
   output.flo : inpainted flow
   guide.png  : guiding image


-------------------------------------------
LAPLACE-BELTRAMI OPTICAL FLOW INTERPOLATION
-------------------------------------------

We provide a Matlab implementation of the Laplace-Beltrami optical flow interpolation.

 - Test:
 
laplace_beltrami_flow_interpolation('test/input.flo', 'test/mask.png', 'test/guide.png', 'test/output_lb.flo', 0.000001, 3)

 - Usage: 

laplace_beltrami_flow_interpolation(input.flo, mask.png, guide.png, output.flo, lambda, g)

Required parameters:
   lambda     :   anisotropic weight
   g          :   weight selection
   input.flo  : input flow to inpaint
   mask.png   : inpainting mask
   guide.png  : guiding image
   output.flo : inpainted flow

---------------------
COPYRIGHT AND LICENSE
---------------------

Copyright (c) 2019 Lara Raad <lara.raad@upf.edu>,
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
along with this program. If not, see <http://www.gnu.org/licenses/>.
