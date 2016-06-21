function mse = MeanSquareError(x,y)
% This function computes the mean square error 
% between x and y 

mse = sum( sum( (y - x).^2 ) );
mse = mse / prod(size(x));
