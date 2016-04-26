%
% This toolbox contains Matlab files that implement the locally % constructed nonsubsampled shearlet transform developed by 
% Glenn R. Easley, Demetrio Labate, and Wang-Q Lim as described % in Ref. [1,2] below. 
% An alternative implementation of the algorithm is
% is provided in Shearlab-1.1. 
%
% References:
%
% 1. G. Easley, D. Labate and W. Lim, 
%    "Optimally Sparse Image Representations using Shearlets",
%    Proc. 40th Asilomar Conf. on Signals, Systems and 
%    Computers, Monterey (2006)
%
% 2. G. Easley, D. Labate and W. Lim,
%    "Sparse Directional Image Representations using the 
%    Discrete Shearlet Transform",
%    Appl. Comput. Harmon. Anal. 25 pp. 25-46, (2008)
%
%
% The functions used to implement the nonsubsampled Laplacian 
% pyramid filter are: atrousdec, atrousc, atrousfilter, 
% atrousrec, upsample2df, and symext.
% These come from the nsct_toolbox found at
% http://www.mathworks.com/matlabcentral/fileexchange/10049
% and were written by Arthur L. Cunha.
%
% References:
%
% 1. Jianping Zhou, Arthur L. Cunha, and Minh N. Do,
%    "Nonsubsampled contourlet transform: construction and application in enhancement",
%    IEEE International Conference on Image Processing, Genoa, Italy, Sep. 2005.
%
% 2. Arthur L. Cunha, Jianping Zhou, and Minh N. Do,
%    "Nonsubsampled contourlet transform: filter design and application in image denoising,"
%    IEEE International Conference on Image Processing, Genoa, Italy, Sep. 2005.
%
% 3. Arthur L. Cunha, Jianping Zhou, and Minh N. Do, 
%    "Nonsubsampled contourlet transform: Theory, design, and applications", 
%    IEEE Transactions on Image Processing, vol. 15, no. 10, pp. 3089-3101, Oct. 2006. 
%




