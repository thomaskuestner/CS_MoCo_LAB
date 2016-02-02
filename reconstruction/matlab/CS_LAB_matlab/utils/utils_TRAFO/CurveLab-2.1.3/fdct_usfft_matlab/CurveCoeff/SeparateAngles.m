function R = SeparateAngles(X,deep,is_real)
% SeparateAngles: Frequency angular partitioning
%  Usage:
%     R = SeparateAngles(X,deep)
%  Inputs:
%    X     m by m matrix (m = 2*2^(j+1)). 
%          The input is the Fourier transform of an object obtained
%          after windowing with a bandpass Meyer window which
%          isolates frequencies in the range 2^(j-1) <= |k| <=
%          2^j. We interpret the array as samples at the frequency
%          points (k1,k2), -2^(j+1) <= k1,k2 < 2^(j+1) 
%
%   deep   2^deep is the number of angular wedges in each of the
%          directional  panels NESW (North, East, South, West).
%  Outputs:
%    R      Array of smoothly localized Fourier samples
%           size = 4 * m * 2^j *  2*2^(j+1)/m, m = 2^deep 
%
%     (i)   first index -- cardinal point 1 = W, 2 = E, 3 = N, 4 = S. 
%     (ii)  second index -- indexes orientations according to the
%              organization below; e.g. (1,3) would correspond to
%              the wedge marked with a star.   
%
%                                j 
%             ------------------------------------->
%            |
%            |  1                                4
%            |  
%            |
%            |
%            |  2                                3
%            |
%       i    |
%            | 
%            |  3*                                2
%            |
%            |
%            |  
%            |  4                                1   
%            \/
%
%
%                     West                      East
%
%              North and South are obtained by transposition. 
%
%     (iii) third and fourth indices-- frequency points, regular
%           samples along parallel slanted lines (shape of a parallelogram)
%
% Description 
%   Smoothly localizes the FT near angular wedges  
%   Operates by interpolating the FT along vertical and horizontal 
%   directions
% See Also
%   Adj_SeparateAngles, AtA, AtA_Toeplitz
%
% By Emmanuel Candes, 2003-2004

  if nargin < 3
    is_real = 0;
  end
  
  n  = size(X,1);
  n2 = n/2;
  boxlen = n/2^deep;
  boxcnt = 2^deep;	
  
  R = zeros(4,boxcnt,n2/2,2*boxlen);
  
  [ix,w] = DetailMeyerWindow([n2/4 n2/2],3);
  lx = reverse(-ix);
  
  m      =  0:(boxcnt - 1);
  ym     =  (m-boxcnt/2).*boxlen + boxlen/2;
  slope  =  -ym./n2; 
  
  % recall that fft_mid0 operates along columns
  % Columns: West and East
  
  F = ifft_mid0(X).*sqrt(n); 
  
  for r = 1:(n2/2),   
    k = lx(r); t = n2 + k + 1;
    shift  =  ym + slope.*(k+n2);
    alpha =  -k./n2;
    w = MakeSineWindow(boxlen,alpha);
    y = Evaluate_FT(F(:,t),shift,boxlen,0,w);    
    R(1,:,r,:)  =  reshape(y,2*boxlen,boxcnt).';
    rsym = n2/2+1-r;
    if (is_real)
      R(2,:,rsym,:) = conj(fliplr(squeeze(R(1,:,r,:))));
    else
      t = n2 - k + 1;
      shift = -shift + 1;
      y = Evaluate_FT(F(:,t),shift,boxlen,0,w);
      R(2,:,rsym,:) = reshape(y,2*boxlen,boxcnt).';
    end
  end
  
  % Rows: North and South
  F = ifft_mid0(X.').*sqrt(n); 
  for r = 1:(n2/2),	       
    k = lx(r);  t = n2 + k + 1;
    shift  =  ym + slope.*(k+n2);
    alpha =  -k./n2;
    w = MakeSineWindow(boxlen,alpha);
    y = Evaluate_FT(F(:,t),shift,boxlen,0,w);
    R(3,:,r,:) = reshape(y,2*boxlen,boxcnt).';
    
    rsym = n2/2+1-r;
    if (is_real)
      R(4,:,rsym,:) = conj(fliplr(squeeze(R(3,:,r,:))));
    else
      t = n2 - k + 1;
      shift = -shift + 1;
      y = Evaluate_FT(F(:,t),shift,boxlen,0,w);
      R(4,:,rsym,:) = reshape(y,2*boxlen,boxcnt).';
    end
  end
