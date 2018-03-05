function b = A_fft2( x, mask, P, m, n, bComplex )
% 2D fft
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

N = length(x);
if(nargin < 3), P = 1:N; end

if(~bComplex)
    % real input image
    fct = 1/sqrt(N);
    fx = zeros(m,n);
    fx(P) = x(P);
else
    % complex input image
    fct = 1/sqrt(N/2);
    fx = complex(zeros(m,n),zeros(m,n));
    fx(P(1:end/2)) = sqrt(2)*x(1:N/2) + 1i *sqrt(2)* x(N/2+1:N);
end

fx = fct * fftshift(fft2(ifftshift(fx)));
b = [sqrt(2)*real(fx(mask)); sqrt(2)*imag(fx(mask))];

end

