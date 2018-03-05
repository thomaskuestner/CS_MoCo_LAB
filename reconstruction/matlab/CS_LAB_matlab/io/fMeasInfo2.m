function iLC = fMeasInfo2(sFilename)
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
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Create the involved filenames
    [sPath, sName] = fileparts(sFilename);
    if isempty(sPath)
        sFullName = which(sFilename);
        if isempty(sFullName), error('File ''%s'' not on path!', sFilename); end
        [sPath, sName, sExt] = fileparts(sFullName);
        sFilename = [sPath, filesep, sName, sExt];
    else
        if ~exist(sFilename, 'file'), error('File ''%s'' doesn''t exist!', sFilename); end
    end
    sFilenameMat = [sPath, filesep, sName, '.mat'];
    load (sFilenameMat);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    fprintf(1, '*** Info of file ''%s'' ***\n', sFilename);
    SDrecksMDH.Seq
    fprintf(1, 'Loop counter summary:\n');
    for iI = 1:length(csLCS)
        if(iI == 1)
            fprintf(1, '  %s: [%u..%u]\n', csLCS{iI}, min(iLC(:,iI)), max(iLC(:,iI)));
        else            
            fprintf(1, '  %s: [%u..%u]\n', csLCS{iI}, min(iLC(iLC(:,1) > 20,iI)), max(iLC(iLC(:,1) > 20,iI)));
        end
    end
    
else
    SDir = dir('*.dat');
    for iI = 1:length(SDir)
        fMeasInfo2(SDir(iI).name);
    end
end

