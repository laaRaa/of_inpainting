%Copyright (c) 2019 Lara Raad <lara.raad@upf.edu>
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU Affero General Public License as
%published by the Free Software Foundation, either version 3 of the
%License, or (at your option) any later version.
%
%You should have received a copy of the GNU Affero General Public License
%along with this program. If not, see <http://www.gnu.org/licenses/>. 


function merge_mask(input_2_path,mask_0_path)

input_2 = imread(input_2_path);
mask_0 = imread(mask_0_path);
input_2 = input_2/max(max(input_2));
mask_0 = mean((mask_0>0)*1,3);
mask = (input_2 | mask_0)*255;
imwrite(mask,input_2_path)
