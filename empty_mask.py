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

def empty_mask(input_0_path,input_2_path):
    frame_0 = Image.open(input_0_path)
    mask = np.zeros((frame_0.height,frame_0.width))
    mask = Image.fromarray(np.uint8(mask))
    mask.save(input_2_path)

empty_mask('input_0.png','input_2.png')
