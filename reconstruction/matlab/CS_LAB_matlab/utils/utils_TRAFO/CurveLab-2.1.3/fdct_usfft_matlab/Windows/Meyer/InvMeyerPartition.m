function x = InvMeyerPartition(y, L, deg)
% InvMeyerPartition: Inverse the scale partitioning
%  Usage:
%     y = MeyerPartition(xhat, L, deg);
%  Inputs: 
%    y      Vector of length 2*n; signal separated into disjoint scales 
%    L      Coarsest Scale
%    deg    Degree of the polynomial
%  Outputs:
%    x     Vector of length n
%  Description
%    Performs the inverse of MeyerPartition; i.e., reconstruct
%    an object from its different scale contributions.
%
% By Emmanuel candes, 2003-2004


        n = length(y)/2;
	J = log2(n);
	x = zeros(1,n);

% 
%  Unfold Partition at Coarse Level.
%
        [index, window] = CoarseMeyerWindow(L-1,deg);
	l_index = reverse(n/2 + 1 - index);
	r_index = n/2 + index;
	x(l_index) = y(1:2^L).* reverse(window);
	x(r_index) = y((2^L+1):2^(L+1)).* window;  
	
	
%
%  Loop to Unfold Partition for  j = L - 1, ..., J - 3.
%
	for j = L-1:(J-3),
	  dyadic_points = [2^j 2^(j+1)];
	  [index, window] = DetailMeyerWindow(dyadic_points,deg); 
	  yy = y((2^(j+2)+1):(2^(j+3)));
	  m = length(yy)/2;
	  l_index = reverse(n/2 + 1 - index);
	  x(l_index) = x(l_index) + yy(1:m).* reverse(window);
	  r_index = n/2 + index;
	  x(r_index) = x(r_index) + yy((m+1):(2*m)).* window; 	   
	end

%
%  Finest Subband (for j = J - 2).
%
        
        j = J - 2;
        [index, window] = FineMeyerWindow(j,deg);
	yy = y((2^(j+2)+1):(2^(j+3)));
	m = length(yy)/2;
	l_index = reverse(n/2 + 1 - index);
        x(l_index) = x(l_index) + yy(1:m).* reverse(window);
        r_index = n/2 + index;
	x(r_index) = x(r_index) + yy((m+1):(2*m)).* window; 	  
	
	
	
	