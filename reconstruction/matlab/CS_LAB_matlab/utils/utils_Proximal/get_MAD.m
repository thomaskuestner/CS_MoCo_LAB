function [ threshold ] = get_MAD( x,nCha,TransformSpecifics )
%GET_MAD Summary of this function goes here
%  Used to acquire lambda dependent on fine-scaled coeffs of (atm) dwt.
%  Other transforms could also be used.

% M. Fischer, April 2016

%%
coeff = cell(1,nCha);
threshold = zeros(1:nCha);

if size(x{1,1},3) == 1
    for j=1:nCha
        coeff{1,j} = wavedec2(abs(x{1,j}),4,'db2'); % atm only based on dwt % TransformSpecifics.waveletStages,TransformSpecifics.waveletFilterName
        coeff_fine_scale = size(coeff{1,j},2) - (3*TransformSpecifics.waveS(TransformSpecifics.Stages+1,1)*TransformSpecifics.waveS(TransformSpecifics.Stages+1,2));
        threshold(j) = mad(coeff{1,j}(coeff_fine_scale:end),1);
        if isnan(threshold(j))
            error('threshold value NaN is a bad idea');
        end;
    end;
else
    for j=1:nCha
        coeff{1,j} = wavedec3(abs(x{1,j}),TransformSpecifics.Stages,'db2'); % atm only based on dwt
        coeff_fine_scale = zeros(0,0);
        for k = (7*(coeff{1,j}.level-1)+2):TransformSpecifics.Stages*7+1;
            coeff_fine_scale = [coeff_fine_scale  coeff{1,j}.dec{k}];
        end;
        [size1, size2, size3] = size(coeff_fine_scale);
        threshold(j) = mad(reshape(coeff_fine_scale,size1*size2*size3,1));
        if isnan(threshold(j))
            error('threshold value NaN is a bad idea');
        end;
    end;
end

