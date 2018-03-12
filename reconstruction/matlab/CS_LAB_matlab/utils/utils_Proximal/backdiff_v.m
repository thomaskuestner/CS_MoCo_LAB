function [ y ] = backdiff_v( X,n1,n2 )
%DIFF_V calculate backpropagation vertical difference
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

% X = reshape(x,n1,n2);
% y = [X(1,:) ; diff(X)];
% y = [zeros(1,n2) ; diff(X)];

% X_temp_diff = diff(x);
% y = [zeros(1,n1) ; abs(X_temp_diff)];

% y = [zeros(1,n1) ; X(1:n1-1,:)-X(2:n1,:)];

% zero boundary condition (GAPG): x(0) = x(n) = 0; -> 
% y(n) = x(n-1), y(1) = - x(1)
% x(n) is already zero, see diff_v.

y = [-X(1,:); X(1:n1-1,:) - X(2:n1,:)];
end

