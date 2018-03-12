function dImgLow = fDownSample( dImgHigh, out_size )

out_x_sz = out_size(1); out_y_sz = out_size(2); out_z_sz = out_size(3);
[in_x_sz,in_y_sz, in_z_sz, nCha, nTime] = size( dImgHigh );


% build grid vectors for the up/down sampling
% ============================================
% if the input is even & output is odd-> use floor for all
% if the output is even & input is odd -> use ceil for all
% other cases - don't care
% for downsampling -> the opposite
if (~mod( in_x_sz,2 ) & (out_x_sz>in_x_sz)) | (mod( in_x_sz,2 ) & (out_x_sz<in_x_sz))
    x_output_space  = max(floor((out_x_sz-in_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
    x_input_space   = max(floor((in_x_sz-out_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
else
    x_output_space  = max(ceil((out_x_sz-in_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
    x_input_space   = max(ceil((in_x_sz-out_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
end
if (~mod( in_y_sz,2 ) & (out_y_sz>in_y_sz)) | (mod( in_y_sz,2 ) & (out_y_sz<in_y_sz))
   y_output_space  = max(floor((out_y_sz-in_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
   y_input_space   = max(floor((in_y_sz-out_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
else
   y_output_space  = max(ceil((out_y_sz-in_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
   y_input_space   = max(ceil((in_y_sz-out_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
end
if (~mod( in_z_sz,2 ) & (out_z_sz>in_z_sz)) | (mod( in_z_sz,2 ) & (out_z_sz<in_z_sz))
   z_output_space  = max(floor((out_z_sz-in_z_sz)/2),0) + [1:min(in_z_sz,out_z_sz)];
   z_input_space   = max(floor((in_z_sz-out_z_sz)/2),0) + [1:min(in_z_sz,out_z_sz)];
else
   z_output_space  = max(ceil((out_z_sz-in_z_sz)/2),0) + [1:min(in_z_sz,out_z_sz)];
   z_input_space   = max(ceil((in_z_sz-out_z_sz)/2),0) + [1:min(in_z_sz,out_z_sz)];
end

dImgLow = zeros( out_x_sz, out_y_sz, out_z_sz, size(dImgHigh,4), size(dImgHigh,5) );
dImgLow(x_output_space,y_output_space, z_output_space, :, : ) = dImgHigh(x_input_space,y_input_space, z_input_space, :, :);

end

