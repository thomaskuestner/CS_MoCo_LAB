function x = At_fft2( b, N, mask, P, m, n, bComplex )
% 2D IFFT
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

if( nargin < 4 ), P = 1:N; end

K = length(b);
fx = complex(zeros(m,n),zeros(m,n));
fx(mask) = 1/sqrt(2)*b(1:K/2) + 1i *1/sqrt(2)* b(K/2+1:K);
x = zeros(N,1);

if(bComplex)
    fct = sqrt(N/2); % complex input image
else
    fct = sqrt(N); % real input image
end

tmp = fct*fftshift(ifft2(ifftshift(fx)));

if(~bComplex)
    % real input image
    x(P) = abs(tmp(:));
else
    % complex input image
    x(P) = [1/sqrt(2)*real(tmp(:)); 1/sqrt(2)*imag(tmp(:))];
end

% x(P) = real(tmp(:)); % phantom

end

