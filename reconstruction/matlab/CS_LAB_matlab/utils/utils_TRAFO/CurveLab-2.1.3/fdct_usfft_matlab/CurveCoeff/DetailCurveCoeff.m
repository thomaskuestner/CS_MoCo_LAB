function D = DetailCurveCoeff(X,deep,is_real);
% DetailCurveCoeff: Returns the "tight" frame coefficients of the 
%                   image at the intermediate scale 2^j;
%  Usage:
%    D = DetailCurveCoeff(X,deep,is_real);
%  Inputs:
%    X     m by m matrix, m = 2*2^j. Typically, X is the matrix of
%          coefficients obtained after separating an object into a
%          series of scales. Scale is explicitely specified by the
%          size of X
% 
%    deep  determines the number of angular wedges per directional
%          panel (cardinal point). The total number of wedges
%          (orientations) is 4 * nw, where nw = 2^deep
%
%  Outputs:
%    D    Cell array of curvelet coefficients at intermediate scale. 
%         D{l}(k1,k2) is the coefficient at
%            - Orientation l: indexes wedge angle starting from top left
%            wedge and increasing in a clockwise fashion.  
%            - Position k1, k2: both integers.  The array index is 
%                 effectively a translation parameter. 
% See Also
%   SeparateScales, SeparateAngles, Inv_DetailCurveCoeff
%
% By Emmanuel Candes, 2003-2004

  R = SqueezeAngularFT(SeparateAngles(X,deep,is_real));
  
  for w=1:size(R,2)
    tmp = squeeze(R(2,w,:,:));
    R(2,w,:,:) = tmp([end,1:end-1], [end,1:end-1]);
  end
  for w=1:size(R,2)
    tmp = squeeze(R(4,w,:,:));
    R(4,w,:,:) = tmp([end,1:end-1], [end,1:end-1]);
  end
  
  nn = size(R);
  D = zeros(nn);
  for j = 1:size(R,1),
    for m = 1:size(R,2),
      W = squeeze(R(j,m,:,:));
      W = ifft2_mid0(W)*sqrt(prod(size(W)));
      D(j,m,:,:) = W;
    end
  end
  
  D = WENStoClockwise(D);
    
