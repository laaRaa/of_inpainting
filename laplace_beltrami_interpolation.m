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

function u = laplace_beltrami_interpolation(v,g,kappa,lambda,d)

% input flow components
v1 = v(:,:,1); v1 = v1(:); % x component
v2 = v(:,:,2); v2 = v2(:); % y component
% guide image colour components
if ismatrix(g)
    g = repmat(g,[1 1 3]);
end
gr = g(:,:,1); gr = gr(:); % red channel
gg = g(:,:,2); gg = gg(:); % green channel
gb = g(:,:,3); gb = gb(:); % blue channel

% inainting mask
kappa = kappa/max(max(kappa)); % mask of ones and zeros
[h,w] = size(kappa); % width and height of inputs

% compute weighted Laplacian Lw
B = incidence_matrix_4N(h,w); % incidence matrix of graph hxw

m = size(B,1); % number of edges
[X,Y] = meshgrid(1:w,1:h);
switch d
    case 1
        D = sqrt( (1-lambda)*( (B*gr).*(B*gr) + (B*gg).*(B*gg) + (B*gb).*(B*gb) )/3 +...
    lambda*( (B*X(:)).*(B*X(:)) + (B*Y(:)).*(B*Y(:)) ) );
    case 2
        D = (1-lambda)*sqrt( (B*gr).*(B*gr) + (B*gg).*(B*gg) + (B*gb).*(B*gb) )/3 +...
    lambda*sqrt( (B*X(:)).*(B*X(:)) + (B*Y(:)).*(B*Y(:)) );
    case 3
        D = (1-lambda)*( (B*gr).*(B*gr) + (B*gg).*(B*gg) + (B*gb).*(B*gb) )/3 +...
    lambda*( (B*X(:)).*(B*X(:)) + (B*Y(:)).*(B*Y(:)) );
    otherwise
        disp('distance identifier not valid');
end
% W = exp(-D); W1. we use the following weight since for this case the
% range very small.max(W1(:))/min(W1(:))=~4 and max(W2(:))/min(W2(:))=~10^6
W = 1./D; %W2
% max(W(:))/min(W(:))
W = spdiags(W,0,m,m); % weight matrix
Lw = -B'*(W*B); % weighted Laplacian operator

% sparse matrix containing the mask in the diagonal
kappa = kappa(:);
kappas = spdiags(kappa==1,0,w*h,w*h);
A = speye(w*h) - kappas + kappas*Lw;
% condest(A)
% issymmetric(A)
b1 = (speye(w*h) - kappas)*v1;
b2 = (speye(w*h) - kappas)*v2;
x1 = A\b1;
x2 = A\b2;
u1 = reshape(x1,h,w);
u2 = reshape(x2,h,w);
u = cat(3,u1,u2);
