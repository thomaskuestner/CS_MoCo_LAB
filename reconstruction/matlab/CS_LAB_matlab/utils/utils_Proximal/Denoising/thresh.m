function [ y ] = thresh( x, lam, a, pen )
% y = thresh(x,lam)
% Generalized thresholding in case of non-convex regulzarition with log/atan penalty
% Ankit Parekh, NYU School of Engineering
% ankit.parekh@nyu.edu

y = zeros(size(x));
ind = abs(x)>=lam;

if ~lam
    y = x;
else

    switch pen
        case 'log'
            y(ind) = (abs(x(ind))/2-1/(2*a) + sqrt((abs(x(ind))./2 + 1/(2*a)).^2 - lam/a)).*sign(x(ind));
        case 'atan'
            y = atanT2(x,lam,a);
        case 'l1'
            y = soft(x,lam);
        otherwise
            disp('Please select penalty from the following: log, atan, l1')
    end
end


