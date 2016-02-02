function res = mtimes(a,b)

if isa(a,'Wavelet') == 0
    error('In  A.*B only A can be Wavelet operator');
end

res = b*0;

for n=1:size(b,3)

if a.adjoint
    res(:,:,n) = IWT2_PO(real(b(:,:,n)),a.wavScale,a.qmf) + 1i*IWT2_PO(imag(b(:,:,n)),a.wavScale,a.qmf);
else
    res(:,:,n) = FWT2_PO(real(b(:,:,n)),a.wavScale,a.qmf) + 1i* FWT2_PO(imag(b(:,:,n)),a.wavScale,a.qmf);
end
end


