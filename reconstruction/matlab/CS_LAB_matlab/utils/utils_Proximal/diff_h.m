function [ y ] = diff_h( X,n1,n2 )
%DIFF_V calculate horizontal difference
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

% X = reshape(x,n1,n2);
% y = [X(1,:) ; diff(X)];
% y = [zeros(1,n2) ; diff(X)];

% X_temp_diff = diff(x);
% y = [zeros(1,n1) ; abs(X_temp_diff)];

% y = [X(2:n1,:)-X(1:n1-1,:) ; zeros(1,n2)];  % zeros, since reflexive boundary condition: x(n+1) = x(n) (GAPG).
y = [X(2:n1,:)-X(1:n1-1,:); zeros(1,n2)]; % 256 -> 256
end

