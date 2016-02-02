function X = Adj_DetailCurveCoeff(C,is_real);
% Adj_DetailCurveCoeff: Adjoint of DetailCurveCoeff  
%  Usage:
%    X = Adj_DetailCurveCoeff(C);
%  Inputs:
%    C    matrix of curvelet coefficients at scale 2^j
%  Outputs:
%    X    matrix of Fourier samples; jth dyadic subband
%  See Also
%   DetailCurveCoeff, Adj_SeparateAngles, Adj_Curvelet02Xform
%
% By Emmanuel Candes, 2003-2004

  C = ClockwisetoWENS(C);
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
  
  X = Adj_SeparateAngles(Adj_SqueezeAngularFT(R),deep,is_real);
  
  
  
  
