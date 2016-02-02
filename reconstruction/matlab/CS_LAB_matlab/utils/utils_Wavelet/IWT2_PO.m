function x = IWT2_PO(wc,L,qmf)
% IWT2_PO -- Inverse 2-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    x = IWT2_PO(wc,L,qmf)
%  Inputs
%    wc    2-d wavelet transform [n by n array, n dyadic]
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    x     2-d signal reconstructed from wc
%
%  Description
%    If wc is the result of a forward 2d wavelet transform, with
%    wc = FWT2_PO(x,L,qmf), then x = IWT2_PO(wc,L,qmf) reconstructs x
%    exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FWT2_PO, MakeONFilter
%
	[n,J] = quadlength(wc);
	x = wc; 
	nc = 2^(L+1);
	for jscal=L:J-1, % from coarse to fine
		top = (nc/2+1):nc; bot = 1:(nc/2); all = 1:nc;
		for iy=1:nc,
			x(all,iy) =  UpDyadLo(x(bot,iy)',qmf)'  ...
					   + UpDyadHi(x(top,iy)',qmf)'; 
		end
		for ix=1:nc,
			x(ix,all) = UpDyadLo(x(ix,bot),qmf)  ... 
					  + UpDyadHi(x(ix,top),qmf);
		end
		nc = 2*nc;
	end
	
%
% Copyright (c) 1993. David L. Donoho
%     
    
    
    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
