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

% writes flow .flo in flow .png
% function []=flo_to_png(flow_flo, flow_png)
% flow_flo: flow .flo path
% flow_png: flow .png path

function flo_to_png(flow_in_flo, flow_out_flo, flow_out_png)

flow_in = readFlowFile(flow_in_flo);
mod = sqrt(flow_in(:,:,1).^2 + flow_in(:,:,2).^2);
max_mod = max(max(mod));

flow = readFlowFile(flow_out_flo);
imwrite(flowToColor(flow,max_mod),flow_out_png)
