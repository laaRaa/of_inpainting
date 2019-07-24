function B = grid_incidence_4N(w, h)                     % grid graph WxH

x = sparse(1:w-1, 2:w, 1, w-1, w) - speye(w-1,w);   % path of length W
y = sparse(1:h-1, 2:h, 1, h-1, h) - speye(h-1,h);   % path of length H
B = [ kron(speye(h),x) ; kron(y,speye(w)) ];  % kronecker union - incidence
% matrix for 4 connectivity
end