function y = MeyerPartition(xhat, L, deg)
% MeyerPartition: Partitions a signal into a series of disjoint scales
%  Usage:
%     y = MeyerPartition(xhat, L, deg);
%  Inputs: 
%    xhat   Input signal
%    L      Coarsest Scale
%    deg    Degree of the polynomial
%  Outputs:
%    y     Vector of length 2*n
%  Description
%    Decomposes the Fourier transform of an object into dyadic
%    subbands by smooth windowing with orthonormal Meyer windows
%
% By Emmanuel candes, 2003-2004


        n = length(xhat);
	J = log2(n);
	y = zeros(1,2*n);
% 
%  Compute Partition at Coarse Level.
%

	[index, window] = CoarseMeyerWindow(L-1,deg);
	left_index = reverse(n/2 + 1 - index);
	left = xhat(left_index).* reverse(window);
	right_index = n/2 + index;
	right = xhat(right_index).* window;
	y(1:2^(L+1)) = [left right];
	
%
%  Loop to Get Partition for  j = L - 1, ..., J - 3.
%
	for j = L-1:(J-3),
	  dyadic_points = [2^j 2^(j+1)];
	  [index, window] = DetailMeyerWindow(dyadic_points,deg); 
	  left_index = reverse(n/2 + 1 - index);
	  left = xhat(left_index).* reverse(window);
	  right_index = n/2 + index;
	  right = xhat(right_index).* window; 
	  y((2^(j+2)+1):(2^(j+3))) = [left right]; 
	end

%
%  Finest Subband (for j = J - 2).
%
        
        j = J - 2;
        [index, window] = FineMeyerWindow(j,deg);
	left_index = reverse(n/2 + 1 - index);
	left = xhat(left_index).* reverse(window);
	right_index = n/2 + index;
	right = xhat(right_index).* window; 
	y((2^(j+2)+1):(2^(j+3))) = [left right]; 
	

	
	
	
	
	









