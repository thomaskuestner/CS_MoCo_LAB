function res = times(a,b)

if isa(a,'Wavelet') == 0
    error('In  A.*B only A can be Wavelet operator');
end

if a.adjoint
    res = IWT2_PO(real(b),a.wavScale,a.qmf) + 1i*IWT2_PO(imag(b),a.wavScale,a.qmf);
else
    res = FWT2_PO(real(b),a.wavScale,a.qmf) + 1i* FWT2_PO(imag(b),a.wavScale,a.qmf);
end


