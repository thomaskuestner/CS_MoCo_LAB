function res = mtimes(a,b)
% res = mtimes(FT, x)
%


if a.adjoint
    bb = reshape(b,a.dataSize(1),a.dataSize(2),a.dataSize(3));
    res = zpad(bb.*a.mask,a.imSize(1),a.imSize(2),a.imSize(3));
    res = ifft2c(res);
    res = res.*conj(a.ph);
    switch a.mode
    	case 0
		res = real(res);
   	case 1
		res = real(res);
    end



else
    bb = reshape(b,a.imSize(1),a.imSize(2),a.imSize(3));
    
    switch a.mode
    	case 0
		bb = real(bb);
   	case 1
		bb = real(bb);
    end
    
    bb = bb.*a.ph; % phase correct
    res = fft2c(bb);
    res = crop(res,a.dataSize(1),a.dataSize(2),a.dataSize(3));
    res = res.*a.mask;
end

if size(b,2) == 1
    res = res(:);
end



    
