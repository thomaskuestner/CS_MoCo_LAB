function u = fJointMCCS( dImg, u, pairs )
% dImg  4D image
% u     initial motion field (row: from, column: to)

SOptions.PyrMax = 4; 
SOptions.PyrMin = 0;
SOptions.PreFilt = 0;        % Pre-filtering
SOptions.MedFilt = 1;        % median filtering
SOptions.Algorithm = 'normal'; % fast
SOptions.Nf = 3; % number of filters used in derivatives
SOptions.Interpolation = 'cubicOMOMS'; % interpolation for transformation

% Parameters for optical flow estimation:
FilterSizes =  2.^(SOptions.PyrMax:-1:SOptions.PyrMin);

parfor iI=1:size(pairs,1)
    fprintf('Registration %d/%d\n', iI, size(pairs,1));
    u{iI,1} = CG_MultiScale_LAP3D(dImg(:,:,:,pairs(iI,2)), dImg(:,:,:,pairs(iI,1)), FilterSizes, SOptions.PreFilt, SOptions.MedFilt, u{iI,1});    
end

end

