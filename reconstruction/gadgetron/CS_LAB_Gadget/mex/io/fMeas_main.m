function [ dData, iLC, iEvalInfoMask, drecksMDH, iLCPositions ] = fMeas_main( sFilename, lReadData, measReadPara )
%FMEAS_MAIN main function for extracting Data, Loop counters, evalInfoMask
% and drecksMDH
%
% (c) Christian Wuerslin, Thomas Kuestner
% ---------------------------------------------------------------------

if(nargin < 1)
    [sFile, sPath] = uigetfile({'*.dat', 'ADC meas files (*.dat)'},'Select measurement file', pwd);
    sFilename = [sPath,filesep,sFile];
end

if(nargin < 2)
    lReadData = true;
end

% input read parameter
if(nargin < 3)
    % Set, Seg, Smp
    measReadPara = {[0 inf], [0 64999], [21 inf]}; % for CS_Trufi and CS_Flash, CS_Flash_Retro, CS_EPI
end

[sPath, sName, sExt] = fileparts(sFilename);

if(~exist(fullfile(sPath,[sName,'.mat']),'file'))
    fprintf('Creating look-up table for measdata file\n');
    fMeasCreateLUT_main(sFilename);
end

% load data
if(lReadData)
    [dData, iLC, iEvalInfoMask] = fMeasRead(fullfile(sPath,[sName,sExt]), 'Set', measReadPara{1}, 'Seg', measReadPara{2}, 'Smp', measReadPara{3});
else
    dData = [];
    iLC = [];
    iEvalInfoMask = [];
end
matfile = whos('-file', fullfile(sPath,[sName,'.mat']));
if(ismember('iLCDrecks', {matfile.name}) && ismember('iLCPositions', {matfile.name})) 
    % drecksMDH 1.0   
    [drecksMDH, iLCPositions] = fMeasReadDrecksMDH(fullfile(sPath,[sName,sExt]));

elseif(ismember('SDrecksMDH', {matfile.name}))
    % drecksMDH 2.0
    load(fullfile(sPath,[sName,'.mat']), 'SDrecksMDH');
    drecksMDH = SDrecksMDH;
    iLCPositions = [];
else
    % no drecksMDH at all
    drecksMDH = [];
    iLCPositions = [];
end

end

