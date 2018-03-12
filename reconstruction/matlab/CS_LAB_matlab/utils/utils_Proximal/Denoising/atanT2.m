function x = atanT2(y,T,a)
% x = atanT(y,T,a)
% 
% THRESHOLDING FUNCTION OF USING ARCTANGENT PENALTY FUNCTION:
%   gives the solution of 
%   x = argmin_x f(x) = 0.5*(y-x)^2 + T*phi(x,a);
%   where
%   phi(x,a) =  2./(a*sqrt(3))*(atan((2*a*abs(x)+1)/sqrt(3)) - pi/6)
%
% INPUT
%   y : data (scalar or multidimensional array)
%   T : threshold (scalar or multidimensional array)
%   a : penalty convexity parameter (a>0)
%       if a is too small (less than 1e-10) there is no benifit of using the 
%       non-convex penalty function, and the result is approximatly equal
%       to using soft-thresholding.
%
% OUTPUT
%   x : output of atan thresholding


% if ( a < 1e-10 )
%     
%     x = soft(y, T);         
%     
% else
    b = abs(y);
    ab = a.*b;
    i = find( ab == 1 );
    
    u = a - ab.*a;
    v = u.^3./(27.*a.^6) - (b - T)./(2.*a.^2) + (u.*(ab - 1))./(6.*a.^4); 
    w = (ab - 1)./(3.*a.^2) + u.^2./(9.*a.^4);
    
    f = ((v.^2 - w.^3).^(1/2) - v).^(1/3);
      
    z = f - u./(3.*a.^2) + w./f;

    if isempty(i) ~= 1
        z(i) = ( (b(i) - T )./a.^2 ) .^ (1/3);
    end

    x = abs(z) .* sign(y) .*( (b-T) >= 0);

% end

end

function x = soft(y, T)
% x = soft(y, T)
%
% SOFT THRESHOLDING
% for real or complex data.
%
% INPUT
%   y : data (scalar or multidimensional array)
%   T : threshold (scalar or multidimensional array)
%
% OUTPUT
%   x : output of soft thresholding
%
% If x and T are both multidimensional, then they must be of the same size.

x = max(1 - T./abs(y), 0) .* y;
end