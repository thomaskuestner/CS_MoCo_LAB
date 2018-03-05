function [ y ] = diff_v( X,n1,n2 )
%DIFF_V calculate vertical difference
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

% X = reshape(x,n1,n2);
% % y = [X(:,1) diff(X,1,2)];
% y = [zeros(n1,1) diff(X,1,2)];

y = [X(:,2:n2) - X(:,1:n2-1) zeros(n1,1)]; % zeros, since reflexive boundary condition: x(n+1) = x(n).

end

