function [drecksMDH, iLCPositions] = fMeasReadDrecksMDH(sFilename)
%FMEASREADDRECKSMDH Read the exported parameters from the Drecks MDH 
% Important: DrecksMDH structure and input (including scaling to fit into 
% unsigned short variable) need to be known! 
%
% Copyright 2013/2014 Thomas Kuestner, University of Tuebingen, Germany
% thomas.kuestner@med.uni-tuebingen.de
%
% here it is:

%%%%%%%%%%%%%%%
%%% 1st MDH %%%
%%%%%%%%%%%%%%%
%%% unsigned short %%%
% 01.)  dimension (3D := 1, 2D := 0)
% 2D:
% 02.)  slices
% 03.)  multislice-mode
% 04.)  multislice-mode interleaving value
% 05.)  slice thickness [mm] * 100
% 06.)  slice distance [mm] * 100
% 3D:
% 02.)  partitions
% 03.)  -
% 04.)  -
% 05.)  slice thickness [mm] * 100
% 06.)  -
% 
% 07.)  Lines (y)
% 08.)  Base resolution (x)
% 09.)  sagittal position [mm] * 10 + 20000 (LR axis)
% 10.)  coronal position [mm] * 10 + 20000 (AP axis)
% 11.)  transversal position [mm] * 10 + 20000 (HF axis)
% 12.)  FOV read [mm]
% 13.)  FOV phase [mm]
% 14.)  phase resolution (y) [%] * 100
%%% float %%%
% 15.)  sagittal norm vector
%%% unsigned short %%%
% 16.)  slice resolution (z) [%] * 100
% 17.)  images per slab
% 18.)  TR [us] => [ms]
% 19.)  TE [us] => [ms]
% 20.)  flip angle
% 21.)  bandwidth per pixel [Hz/px] -> Attention: Siemens UI value is rounded to next decadic position
% 22.)  RF duration [us]
% 23.)  RF Time*Bandwidth * 10
% 24.)  CS Acceleration * 10
% 25.)  CS Fully Sampled [%] * 10
%%%%%%%%%%%%%%%
%%% 2nd MDH %%%
%%%%%%%%%%%%%%%
% 01.)  CS acquired lines
% 02.)  CS fixed lines
% 03.)  CS sampling type
% 04.)  CS VD map
% 05.)  CS additional parameter 1 * 100
% 06.)  CS additional parameter 2
% 07.)  averages
% 08.)  concatenations
% 09.)  slice oversampling * 1000
% 10.)  phase oversampling * 1000
% 11.)  readout oversampling
% 12.)  ESPReSSo direction
% 13.)  ESPReSSo factor
% 14.)  Image lines
%%% float %%%
% 15.)  coronal norm vector
%%% unsigned short %%%
% 16.)  Image columns
% 17.)  partitions (kSpace samples in z direction without correction of slice pf)
% 18.)  partition fft length
% 19.)  lines (kSpace samples in y direction without correction of phase pf)
% 20.)  phase fft length
% 21.)  echo line
% 22.)  echo partition
% 23.)  sequence type (0: CS_FLASH, 21: CS_Trufi 42: CS_EPI)
% 24.)  temporal phases
%%%%%%%%%%%%%%%
%%% 3rd MDH %%%
%%%%%%%%%%%%%%%
%%% float %%%
% 15.)  transversal norm vector

% additional parameters extracted from MDH
% 01.) slice data position vector (sagittal)
% 02.) slice data position vector (coronal)
% 03.) slice data position vector (transversal)
% 04.) patient table position [mm] + 5 ) * -10

%   2013 Thomas Kuestner, University of Tuebingen, Germany
%   thomas.kuestner@med.uni-tuebingen.de

if nargin == 0, error('No filename given.'); end

% Load the loop counters
[sPath, sName] = fileparts(sFilename);
if isempty(sPath)
    load([sName, '.mat']);
else
    load([sPath, filesep, sName, '.mat']);
end

% get amount of channels for stepsize
nCha = max(unique(double(iLC(:,2))));
lmask = iLC(:,1) < 20;

if(exist('iLCPositions','var') && ~isempty(iLCPositions))
    iLCPositions = iLCPositions(~lmask,:);
    iLCPositions = unique(iLCPositions,'rows'); % slice data position vector (sagittal, coronal, transversal shift [mm] out of isocenter and table position [mm])
else
    iLCPositions = [];
end
clear 'iLC' 'iSP';
if(exist('iLCDrecks','var') && ~isempty(iLCDrecks))
    drecksMDH = zeros(size(iLCDrecks,1)/nCha,25,'uint16');
    for i=1:size(drecksMDH,1)
        drecksMDH(i,:) = iLCDrecks((i-1)*nCha+1,[3:16,21,23:32]);
    end
else
    drecksMDH = [];
end

end



































