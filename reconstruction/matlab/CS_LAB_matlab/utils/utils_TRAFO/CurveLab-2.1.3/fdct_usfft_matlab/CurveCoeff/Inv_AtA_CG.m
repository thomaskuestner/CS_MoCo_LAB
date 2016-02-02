function Y = Inv_AtA_CG(X,deep,MaxIts,ErrTol,X0)
% Inv_SeparateAngles: Invert AtA
%  Usage:
%    Y = Inv_AtA_CG(X,deep,MaxIts,ErrTol);
%  Inputs:
%    X       m by m array of Fourier coefficients
%    deep    Depth of the angular splitting
%    MaxIts  Maximum number of CG iterations. Default 10.
%    ErrTol  Error tolerance. Default 1.e-9.
%  Outputs:
%    X   m by m matrix  
%  Description 
%    Performs the inverse AtA. This is an approximate
%    inverse which uses a conjugate gradient solver. Each iteration
%    uses AtA_Toeplitz which is rapidly evaluated as it merely
%    involves 2 FFTs.
% See Also
%   Inv_SeparateAngles, AtA_Toeplitz 
%
% By Emmanuel Candes, 2003-2004, developped after Amir Avebuch's CG
% solver.  
  
  if nargin < 4,
    ErrTol = 1.e-9;
  end
  if nargin < 3,
    MaxIts = 10;
  end
  
  % Initialization
  xk = X;
  
  n = size(X,1);
  Lambda = MakeFourierDiagonal_2D(n,deep);
  
  S = SetScaleToZero(n);
  
  if MaxIts>1,
    temp = AtA_Toeplitz(xk,Lambda).*S;
    
    pk=xk-temp;
    rk=xk-temp;
    
    m=length(xk);
    
    % Conjugate gradient iteration
    for j=2:MaxIts,
      perr=max(max(abs(pk)));
 %     fprintf(['Iteration ',int2str(j-1),': Residual error=']);
 %     fprintf('%f\n',norm(pk(:)));
      
      if perr>ErrTol,			
        temp = AtA_Toeplitz(pk,Lambda).*S;
        a0 = sum(abs(rk(:)).^2);
        a1 = sum(conj(pk(:)).*temp(:));
        
        a=a0/a1;
        xkp1=xk+a*pk;		
        rkp1=rk-a*temp;
        
        bb = sum(abs(rkp1(:).^2));
        b=bb/a0;
        
        % Update
        pk=rkp1+b*pk;
        rk=rkp1;
        xk=xkp1;
        
        if prod(size(X0)) ~= 0,	  
          rel = norm(xk.*S - X0,'fro')/norm(X0,'fro')*1e3;
   %       fprintf(['Iteration ',int2str(j-1),': Relative error=']);
   %       fprintf('%f\n',rel);
        end
      end
    end
  end
  
  Y = xk;
  
