% Copyright (c) 2019 Lara Raad <lara.raad@upf.edu>,
%			       <enric.meinhardt@cmla.ens-cachan.fr>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
%
% Smooth Contours is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

% laplace beltrami flow interpolation
% function []=laplace_beltrami_flow_interpolation(input, output, lambda, g)
% input_flow: texture sample path
% mask: inpainting mask path
% guide: guiding image path
% output_flow: synthesized texture path
% lambda: anisotropic weight
% g: metric

function laplace_beltrami_flow_interpolation(flow_in, mask, guide,...
    flow_out, lambda, g)

kappa = im2double(imread(mask)); % read inpainting mask
kappa = kappa(:,:,1); % 
kappa(kappa~=0) = 1; % values ~= 0 are pixels to inpaint
I = im2double(imread(guide)); % read guiding image
v = readFlowFile(flow_in); % read ground truth flow
u = inpaint_flow(v,I,kappa,lambda,g); % inpaint missing regions
writeFlowFile(u,flow_out)
