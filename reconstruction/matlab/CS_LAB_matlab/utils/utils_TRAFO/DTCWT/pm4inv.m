function [u, v, q, r] = pm4inv(a, b, c, d)

% [u, v, q, r] = pm4(a, b, c, d);
% u = (a + b - c - d)/2;
% v = (a - b + c + d)/2; 
% q = (a + b - c + d)/2;
% r = (a + b + c - d)/2;

u = ( a + b + c + d)/2;
v = (-a - b + c + d)/2; 
q = (-a + b - c + d)/2;
r = (-a + b + c - d)/2;
