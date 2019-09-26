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