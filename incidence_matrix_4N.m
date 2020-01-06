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

function B = incidence_matrix_4N(h, w)                     % grid graph hxw

x = sparse(1:h-1, 2:h, 1, h-1, h) - speye(h-1,h);   % path of length h
y = sparse(1:w-1, 2:w, 1, w-1, w) - speye(w-1,w);   % path of length w
B = [ kron(speye(w),x) ; kron(y,speye(h)) ];  % kronecker union - incidence
% matrix for 4 connectivity
end
