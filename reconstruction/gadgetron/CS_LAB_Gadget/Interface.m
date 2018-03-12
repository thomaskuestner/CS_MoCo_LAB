%% Interface.m
% small example
%   - load k-space
%   - load parameter file
%   - perform reconstruction with Gadgetron based implementation
%
% (c) 2015 Martin Schwartz and Thomas Kuestner
% -------------------------------------------------------------------------

addpath(genpath(pwd));
    
% load k-space from raw Siemens measurement file
[dKSpace, measPara, file] = fLoadkSpace(PATH_TO_SIEMENS_MEAS);
% OR: 
% load k-space from file/workspace
% required format:
% 2D: xFreq x yPhase x (nTime) x nChannels
% 3D: xFreq x yPhase x zPhase x nChannels
% 4D: xFreq x yPhase x zPhase x nTime x nChannels (not yet supported)

% load config XML parameters
sTrafo = fLoadXML('CS_LAB.xml', measPara);

% process data in Gadgetron implementation
dImage = CS_LAB_Interface(dKSpace, sTrafo);