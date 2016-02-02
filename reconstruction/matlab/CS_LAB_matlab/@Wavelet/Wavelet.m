function res = Wavelet(filterType, filterSize, wavScale)
% res = Wavelet(Filtertype, filterSize, wavScale)
%
% implements a wavelet operator. This is a wrapper to David Donoho's 
% Wavelab. 
%
% Inputs:
%		filterType:   string, 'Haar', 'Beylkin', 'Coiflet', 'Daubechies',
%            			'Symmlet', 'Vaidyanathan','Battle'
%                       Use suffix _TI for translation invariance, for
%                       example, 'Daubechies_TI'
%		filterSize: related to the support and vanishing moments of the particular
%				wavelet (See MakeONFilter in wavelab)
%		wavScale: 	smallest scale of wavelet decomposition
%       TIflag: 0 orthogonal, 1 translation invariant
%
%
%
%
% (c) Michael Lustig 2007

[filterType,TI]=strtok(filterType,'_');
if isempty(TI)==0
    res.TI = 1;
else
    res.TI = 0;
end

res.adjoint = 0;
res.qmf = MakeONFilter(filterType, filterSize);
res.wavScale = wavScale;
res = class(res,'Wavelet');
