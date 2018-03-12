function [ drecksMDH, iLCPositions ] = fMeas_mainLUTDrecks( sFilename )
%FMEAS_MAIN main function for extracting Data, Loop counters, evalInfoMask
% and drecksMDH
%
% (c) Christian Wuerslin, Thomas Kuestner
% ---------------------------------------------------------------------

% input read parameter
% Set, Seg, Smp
% measReadPara = {[0 inf], [0 inf], [21 inf]}; % for CS_Trufi and CS_Flash, CS_Flash_Retro, CS_EPI

[sPath, sName, sExt] = fileparts(sFilename);

if(~exist(fullfile(sPath,[sName,'.mat']),'file'))
    fprintf('Creating look-up table for measdata file\n');
    fMeasCreateLUT_main(sFilename);
end

% load data
matfile = whos('-file', fullfile(sPath,[sName,'.mat']));
if(ismember('iLCDrecks', {matfile.name}) && ismember('iLCPositions', {matfile.name})) 
    % drecksMDH 1.0
%     [dData, iLC, iEvalInfoMask] = fMeasRead(fullfile(sPath,[sName,sExt]), 'Set', measReadPara{1}, 'Seg', measReadPara{2}, 'Smp', measReadPara{3});
    [drecksMDH, iLCPositions] = fMeasReadDrecksMDH(fullfile(sPath,[sName,sExt]));

elseif(ismember('SDrecksMDH', {matfile.name}))
    % drecksMDH 2.0
%     [dData, iLC, iEvalInfoMask] = fMeasRead(fullfile(sPath,[sName,sExt]), 'Set', measReadPara{1}, 'Seg', measReadPara{2}, 'Smp', measReadPara{3});
    load(fullfile(sPath,[sName,'.mat']), 'SDrecksMDH');
    drecksMDH = SDrecksMDH;
    iLCPositions = [];
else
    % no drecksMDH at all
%     [dData, iLC, iEvalInfoMask] = fMeasRead(fullfile(sPath,[sName,sExt]), 'Set', measReadPara{1}, 'Seg', measReadPara{2}, 'Smp', measReadPara{3});
    drecksMDH = [];
    iLCPositions = [];
end

end

