function b = A_fft3( x, mask, P, m, n, o, bComplex )
% 3D fft
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

N = length(x);
if(nargin < 3), P = 1:N; end

if(~bComplex)
    % real input image
    fct = 1/sqrt(N);
    fx = zeros(m,n,o);
    fx(P) = x(P);
else
    % complex input image
    fct = 1/sqrt(N/2);
    fx = complex(zeros(m,n,o),zeros(m,n,o));
    fx(P(1:end/2)) = sqrt(2)*x(1:N/2) + 1i *sqrt(2)* x(N/2+1:N);
end

for i=1:3
    fx = fftshift(fft(ifftshift(fx,i),[],i),i);
end
b = fct .* [sqrt(2)*real(fx(mask)); sqrt(2)*imag(fx(mask))];

end

