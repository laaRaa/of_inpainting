#!/usr/bin/env python3

import os, glob, math
import numpy as np
from PIL import Image

#%Copyright (c) 2019 Lara Raad <lara.raad@upf.edu>
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.
#
#You should have received a copy of the GNU Affero General Public License
#along with this program. If not, see <http://www.gnu.org/licenses/>. 

def merge_mask(input_2_path,mask_d_path,mask_path):
    mask_d = Image.open(mask_d_path) 
    mask_d = np.array(mask_d)
    if (mask_d.ndim > 2):
        mask_d = np.mean(mask_d,2)
    input_2 = Image.open(input_2_path)
    input_2 = np.array(input_2)
    input_2 = input_2/np.maximum(1,np.amax(input_2))
    mask = np.logical_or(input_2,mask_d)*255
    mask = Image.fromarray(np.uint8(mask))
    mask.save(mask_path)

merge_mask('input_2.png','mask_d.png','mask.png')
