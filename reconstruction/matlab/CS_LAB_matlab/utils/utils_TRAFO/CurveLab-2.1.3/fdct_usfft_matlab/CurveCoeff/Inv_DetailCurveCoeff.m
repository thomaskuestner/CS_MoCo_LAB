function X = Inv_DetailCurveCoeff(C,IsImageReal);
% Inv_DetailCurveCoeff:  Reconstruct the jth dyadic frequency band
%                       from curvelet coefficients at that scale
%  Usage:
%    X = Inv_DetailCurveCoeff(C);
%  Inputs:
%    C    matrix of curvelet coefficients at scale 2^j
%  Outputs:
%    X    jth dyadic subband
%  See Also
%   DetailCurveCoeff, Inv_SeparateAngles, Inv_Curvelet02Xform
%
% By Emmanuel Candes, 2003-2004

  C =  ClockwisetoWENS(C);
  nn = size(C); 
  R = zeros(nn);
  deep = log2(nn(2));
  
  for j = 1:size(R,1),
    for m = 1:size(R,2),
      W = squeeze(C(j,m,:,:));
      W = fft2_mid0(W)/sqrt(prod(size(W)));
      R(j,m,:,:) = W;
    end
  end
  
  for w=1:size(R,2)
    tmp = squeeze(R(2,w,:,:));
    R(2,w,:,:) = tmp([2:end,1], [2:end,1]);
  end
  for w=1:size(R,2)
    tmp = squeeze(R(4,w,:,:));
    R(4,w,:,:) = tmp([2:end,1], [2:end,1]);
  end
  
  MaxIts = 25; 
  X = Inv_SeparateAngles(Adj_SqueezeAngularFT(R),deep, ...
                         IsImageReal,MaxIts,1e-6,[]); 
