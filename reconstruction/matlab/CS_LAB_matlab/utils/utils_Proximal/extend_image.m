function [ z_extend ] = extend_image( z, wavecoeffS_proxA, waveletStages, proxA_extend_y, proxA_extend_x )
% mirrored extension by proxA_extend_y / proxA_extend_x
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

z_extend = zeros(wavecoeffS_proxA(waveletStages+2,1),wavecoeffS_proxA(waveletStages+2,2));
z_extend(1:end-proxA_extend_y,1:end-proxA_extend_x) = z;
z_extend(end-proxA_extend_y+1:end,1:end-proxA_extend_x) = flipud(z(end-proxA_extend_y+1:end,:));
z_extend(1:end-proxA_extend_y,end-proxA_extend_x+1:end) = fliplr(z(:,end-proxA_extend_x+1:end));
z_extend(end-proxA_extend_y+1:end,end-proxA_extend_x+1:end) = rot90(z(end-proxA_extend_y+1:end,end-proxA_extend_x+1:end),2);

end

