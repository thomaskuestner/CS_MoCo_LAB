function nu = WindowMeyer(xi,deg)
% WindowMeyer -- auxiliary window function for Meyer wavelets.
%  Usage
%    nu = WindowMeyer(xi,deg)
%  Inputs
%    xi     abscissa values for window evaluation
%    deg    degree of the polynomial defining Nu on [0,1]
%           1 <= deg <= 3
%  Outputs
%    nu     polynomial of degree 'deg' if x in [0,1]
%           1 if x > 1 and 0 if x < 0.
%  See Also
%    MeyerPartition

	if deg == 0,
		nu = xi;
	else
	  if deg == 1,
		nu = xi .^2 .* (3 - 2 .*xi) ;
	  else
		 if deg == 2,
			nu = xi .^3 .* (10 - 15 .* xi + 6 .* xi .^2);
		 else
			if deg == 3,
				 nu = xi.^4 .* ( 35 - 84 .* xi + 70 .* xi.^2 - 20 .* xi.^3);
			 end
		end
	  end
	end
	ix0 = find(xi <= 0);
	if length(ix0) > 0,
		%size(ix0),
		nu(ix0) = zeros(1,length(ix0));
	end
	ix1 = find(xi >= 1);
	if length(ix1) > 0,
		nu(ix1) = ones(1,length(ix1));
	end


%  Prepared for the thesis of Eric Kolaczyk, Stanford University, 1994
%  Copyright (c) Eric Kolaczyk, 1994.

% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
   
    
