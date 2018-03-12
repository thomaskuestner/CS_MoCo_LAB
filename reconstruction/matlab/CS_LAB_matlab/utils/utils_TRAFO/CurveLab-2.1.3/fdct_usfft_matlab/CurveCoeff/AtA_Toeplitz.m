function Y = AtA_Toeplitz(X,Lambda);
% AtA_Toeplitz: Gram Operator of SeparateAngle
%  Usage:
%     Y = AtA(X,Lambda)
%  Inputs:
%    X      n*n matrix (x,y)
%    Lambda Fourier multipliers
%  Outputs:
%    Y      n*n matrix (x,y)
%  Description:
%    Performs  A'A where A is the angular windowing using the
%    Toeplitz structure of A'A. AtA_Toeplitz uses FFTs.
%
% By Emmanuel Candes, 2003-2004

  n  = size(X,1); n2 = n/2;
  
  X1 = zeros(n); X2 = zeros(n);
  
  [ix,w] = DetailMeyerWindow([n2/4 n2/2],3);
  lx = reverse(-ix);
  
  F = ifft_mid0(X).*sqrt(n); % Take FT
  
  for r = 1:(n2/2),
    k = lx(r); t = n2 + k + 1;
    X1(:,t) = GUSFT_Toeplitz(F(:,t),Lambda(:,r,1)); 
    
    rsym = n2/2+1-r; t = n2 - k + 1;
    X1(:,t) = GUSFT_Toeplitz(F(:,t),Lambda(:,rsym,2)); 
  end
  
  X1 = fft_mid0(X1)/sqrt(n); % Invert FT
  
  F = ifft_mid0(X.').*sqrt(n); % Take FT 
  
  for r = 1:(n2/2),
    k = lx(r); t = n2 + k + 1;
    X2(:,t) = GUSFT_Toeplitz(F(:,t),Lambda(:,r,1)); 
    
    rsym = n2/2+1-r;  t = n2 - k + 1;
    X2(:,t) = GUSFT_Toeplitz(F(:,t),Lambda(:,rsym,2)); 
  end
  
  X2 = fft_mid0(X2)/sqrt(n); % Invert FT
  
  Y = X1 + X2.';
  
  
