function kSpaceLow = updownsample( kSpaceHigh, lInputImage, out_size, lAAfilter)
%
% updownsample - up-sample or down-sample an input series using fourier domain
%                input series needs to be continuous of a high degree
%
%
% input:    in_m                - input matrix for up/down sampling, in fourier domain
%                                 
%           out_size            - desired number of pixels in the output image

%
% output:   out_m               - up/down sampled image 
%
% NOTE: it is important to specify if the image is REAL or COMPLEX, since
%       if the image is REAL -> we have to use ABS() on the inverse fourier
%       transform (because of roundoff errors of the transform).
%
% NOTE: since a desired amount of pixels is needed at the output, there is
%       no attempt to use matrices which are in size of power of 2. this
%       optimization can not be used in this case
%
% NOTE: input series needs to be CONTINUOUS OF A HIGH DEGREE, since the
%       upsampling is done in the frequency domain, which samples the output
%       grid with SINE-like (harmonic) functions
%
%
% 
% Theory:   the upsampling is done by zero-padding in the input domain BETWEEN the samples,
%           then at the fourier domain, taking the single spectrum (out of the repetition of spectrums)
%           i.e. low pass with Fcutoff=PI/upsample_factor, zeroing the rest of the spectrum
%           and then doing ifft to the distribution.
%           since we have a zero padding operation in time, we need to multiply by the fourier gain.
%
%              +-----------+     +-------+     +---------+     +--------+     +--------+
%   y[n,m] --> | up-sample | --> |  FFT  | --> |   LPF   | --> | * Gain | --> |  IFFT  | --> interpolated 
%              | factor M  |     +-------+     | Fc=PI/M |     +--------+     +--------+
%              +-----------+                   +---------+
%
%           this operation is the same as the following one (which has less operations):
%
%              +-------+     +--------+     +--------------+     +--------+
%   y[n,m] --> |  FFT  | --> | * Gain | --> | Zero Padding | --> |  IFFT  | --> interpolated 
%              +-------+     +--------+     +--------------+     +--------+
%
%           NOTE THAT, the zero-padding must be such that the D.C. ferquency remains the D.C. frequency
%           and that the zero padding is applied to both positive and negative frequencies. 
%           The zero padding actually condences the frequency -> which yields a longer series in the 
%           image domain, but without any additional information, thus the operation must be an interpolation.

if(nargin < 4)
    lAAfilter = true;
end

% ==============================================
% get input image size, and calculate the gain
% ==============================================
out_x_sz = out_size(1); out_y_sz = out_size(2); out_z_sz = out_size(3);
[in_x_sz,in_y_sz, in_z_sz, nCha, nTime] = size( kSpaceHigh );
gain_x = out_x_sz/in_x_sz;
gain_y = out_y_sz/in_y_sz;
gain_z = out_z_sz/in_z_sz;

% upsample or downsample as needed
% ==================================
if(lInputImage)
    % image -> kSpace
    kSpaceHigh = fftnshift(kSpaceHigh,1:3);
end

% check if up/down sampling is needed at all
if(gain_x == 1 && gain_y == 1 && gain_z == 1)
    kSpaceLow = kSpaceHigh;
    dImgLow = ifftnshift(kSpaceLow,1:3);
    return;
end

% if downsampling -> apply anti-aliasing filter
if(gain_x < 1 && gain_y < 1 && gain_z < 1 && lAAfilter)
    kSpaceHigh = kSpaceHigh .* repmat(windowND(@hamming, [size(kSpaceHigh,1), size(kSpaceHigh,2), size(kSpaceHigh,3)]),[1 1 1 size(kSpaceHigh,4) size(kSpaceHigh,5)]);
end

% build grid vectors for the up/down sampling
% ============================================
% if the input is even & output is odd-> use floor for all
% if the output is even & input is odd -> use ceil for all
% other cases - don't care
% for downsampling -> the opposite
if (~mod( in_x_sz,2 ) && (out_x_sz>in_x_sz)) || (mod( in_x_sz,2 ) && (out_x_sz<in_x_sz))
    x_output_space  = max(floor((out_x_sz-in_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
    x_input_space   = max(floor((in_x_sz-out_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
else
    x_output_space  = max(ceil((out_x_sz-in_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
    x_input_space   = max(ceil((in_x_sz-out_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
end
if (~mod( in_y_sz,2 ) && (out_y_sz>in_y_sz)) || (mod( in_y_sz,2 ) && (out_y_sz<in_y_sz))
   y_output_space  = max(floor((out_y_sz-in_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
   y_input_space   = max(floor((in_y_sz-out_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
else
   y_output_space  = max(ceil((out_y_sz-in_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
   y_input_space   = max(ceil((in_y_sz-out_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
end
if (~mod( in_z_sz,2 ) && (out_z_sz>in_z_sz)) || (mod( in_z_sz,2 ) && (out_z_sz<in_z_sz))
   z_output_space  = max(floor((out_z_sz-in_z_sz)/2),0) + [1:min(in_z_sz,out_z_sz)];
   z_input_space   = max(floor((in_z_sz-out_z_sz)/2),0) + [1:min(in_z_sz,out_z_sz)];
else
   z_output_space  = max(ceil((out_z_sz-in_z_sz)/2),0) + [1:min(in_z_sz,out_z_sz)];
   z_input_space   = max(ceil((in_z_sz-out_z_sz)/2),0) + [1:min(in_z_sz,out_z_sz)];
end


% perform the up/down sampling
kSpaceLow  = zeros( out_x_sz, out_y_sz, out_z_sz, size(kSpaceHigh,4), size(kSpaceHigh,5) );
% kSpaceHigh = fftshift(fftshift(fftshift(kSpaceHigh,1),2),3); % already zero-centered
kSpaceLow(x_output_space,y_output_space, z_output_space, :, :) = kSpaceHigh(x_input_space,y_input_space, z_input_space, :, :);
    
        
end