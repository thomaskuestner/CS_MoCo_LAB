function [af, sf] = farras

% Farras nearly symmetric filters for orthogonal
% 2-channel perfect reconstruction filter bank
%
% USAGE:
%    [af, sf] = farras
% OUTPUT:
%    af - analysis filters
%    sf - synthesis filters
% REFERENCE:
%    A. F. Abdelnour and I. W. Selesnick. 
%    "Nearly symmetric orthogonal wavelet bases",
%    Proc. IEEE Int. Conf. Acoust., Speech,
%    Signal Processing (ICASSP), May 2001.
% See afb, dwt.
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

af = [
                  0  -0.01122679215254
                  0   0.01122679215254
  -0.08838834764832   0.08838834764832
   0.08838834764832   0.08838834764832
   0.69587998903400  -0.69587998903400
   0.69587998903400   0.69587998903400
   0.08838834764832  -0.08838834764832
  -0.08838834764832  -0.08838834764832
   0.01122679215254                  0
   0.01122679215254                  0
];
 
sf = af(end:-1:1, :);

