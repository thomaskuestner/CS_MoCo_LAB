function iLC = fMeasInfo(sFilename)
% FMEASINFO(SFILENAME) prints loop counter information of the given
% measdata file.
%   
%   FMEASINFO looks for '.dat'-files in the current working directory
%   and prints loop counter information for these files provided the LUT
%   was created using FMEASCREATELUT.
%
%   Class support for input sFilename:
%      string
%
%   See also FMEASCREATELUT, FMEASREAD
%
%   Copyright 2012 Christian Wuerslin, University of Tuebingen, Germany
%   christian.wuerslin@med.uni-tuebingen.de
%   Modifications 2013 Thomas Kuestner, University of Tuebingen, Germany
%   thomas.kuestner@med.uni-tuebingen.de

csLCS = {'Smp', 'Cha', 'Lin', 'Acq', 'Sli', 'Par', 'Eco', 'Pha', 'Rep', 'Set', 'Seg', 'FreePara1', 'FreePara2', 'FreePara3', 'FreePara4'};

if nargin
    [sPath, sName] = fileparts(sFilename);
    if isempty(sPath)
        load([sName, '.mat']);
    else
        load([sPath, filesep, sName, '.mat']);
    end
    
    fprintf(1, 'Loop counter summary:\n');
    for iI = 1:length(csLCS)
        fprintf(1, '  %s: [%u..%u]\n', csLCS{iI}, min(iLC(:,iI)), max(iLC(:,iI)));
    end
    
else
    SDir = dir('*.dat');
    for iI = 1:length(SDir)
        fMeasInfo(SDir(iI).name);
    end
end

