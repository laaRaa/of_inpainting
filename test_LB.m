% Copyright (c) 2019 Lara Raad <lara.raad@upf.edu>,
% Enric Meinhardt-Llopis <enric.meinhardt@cmla.ens-cachan.fr>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

inpaint_OF_with_laplace_beltrami_flow_interpolation('test/input.flo',...
    'test/guide.png', 'test/mask.png', 0.000001, 3, 'test/output_lb.flo')
