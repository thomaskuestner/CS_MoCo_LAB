function fMeasCreateLUT_main(sFilename)
%FMEASCREATELUT Creates look-up table from measdata file(s).
%   FMEASCREATELUT(SFILENAME) parses the given Siemens measdata
%   file and creates a look-up table of the file's MDH entries and saves it
%   to a file '%sFilename%.mat'. This information can be used to load
%   measdata from that file later on using fMeasRead way faster than
%   usually. Additionally, only MDHs with specific loop counter values can
%   be read.
%   
%   FMEASCREATELUT looks for '.dat'-files in the current working directory
%   and creates look-up tables/files for all Siemens measdata files.
%
%   Class support for input sFilename:
%      string
%
%   See also FMEASINFO, FMEASREAD
%
%   Copyright 2014 Christian Wuerslin, University of Tuebingen, Germany
%   christian.wuerslin@med.uni-tuebingen.de
%   Modifications 2014 Thomas Kuestner, University of Tuebingen, Germany
%   thomas.kuestner@med.uni-tuebingen.de

% csLCS = {'Smp', 'Cha', 'Lin', 'Acq', 'Sli', 'Par', 'Eco', 'Pha', 'Rep', 'Set', 'Seg', 'FreePara1', 'FreePara2', 'FreePara3', 'FreePara4'};
% Loop Counter names
%
% (c) Christian Wuerslin, Thomas Kuestner
% ---------------------------------------------------------------------

if nargin
    
    % ---------------------------------------------------------------------
    % Filename give => parse file
    
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
    sFilenameXML = [sPath, filesep, sName, '.xml'];
    sFilenameMat = [sPath, filesep, sName, '.mat'];
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Get filesize and extract the header
    SDir = dir(sFilename);
    iFileSize = SDir(1).bytes;
    fidmeas = fopen(sFilename, 'r');
    iHeaderLength = fread(fidmeas, 1, 'uint32');
    fprintf(1, '\n*** Parsing %s ***\n', sFilename);
    fprintf(1, 'File Size            : %u  (%u MB)\n', iFileSize, round(iFileSize/(1024^2)));
    fprintf(1, 'Length of the Header : %u\n', iHeaderLength);
    
    sXML = fread(fidmeas, iHeaderLength, '*char');
    fidXML = fopen(sFilenameXML, 'w');
    fwrite(fidXML, sXML);
    fclose(fidXML);
    fclose(fidmeas);
    fprintf(1, 'Header exported to file ''%s''\n', sFilenameXML);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Get the MDHs
    [dMDH, iCount] = fMeasCreateLUT_mex(sFilename);
    iCount = uint32(iCount);
    dMDH = dMDH(:,1:iCount)';
    fprintf(1, 'Number of ADCs: %u\n', iCount);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Get the DrecksMDH as struct
    lDrecksMDH = dMDH(:, 4) == 20;
    if(nnz(lDrecksMDH) == 0)
        lVer = -1;
    else
        lVer = bitget(uint16(dMDH(find(lDrecksMDH, 1, 'first'), 6)), 16);
    end
    if  lVer == -1 % no DrecksMDH
        iLCDrecks = [];
        iLCPositions = [];
    elseif lVer == 0 % or dMDH(1,6) ?, bitget(uint16(dMDH(1, 1)), 16) == 0
        % DrecksMDH 1.0
        iLCDrecks = [dMDH(lDrecksMDH,4:19),zeros(nnz(lDrecksMDH),4),dMDH(lDrecksMDH,30),zeros(nnz(lDrecksMDH),1),dMDH(lDrecksMDH,20:29),dMDH(lDrecksMDH,31:37),zeros(nnz(lDrecksMDH),1),dMDH(lDrecksMDH,38)];
        iLCPositions = dMDH(~lDrecksMDH, 31:end);
    else
        % DrecksMDH 2.0
        dDrecksMDH = dMDH(lDrecksMDH, :);
        iNChannels = dDrecksMDH(1, 5);
        dDrecksMDH = dDrecksMDH(1:iNChannels:end, :);
        dInt = dDrecksMDH(:, 6:28)';
        dFloat = dDrecksMDH(:, 29:30)';
        dDrecksMDH = [dInt(:); dFloat(:)];
        SDrecksMDH = fReadNewDrecksMDH(dDrecksMDH);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    iSP = uint64(dMDH(~lDrecksMDH,1));
    iEvalInfoMask = uint32(dMDH(~lDrecksMDH, 2:3));
    iLC = uint16(dMDH(~lDrecksMDH, [4:14,26:29]));
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (lVer == -1)
        save(sFilenameMat, 'iLC', 'iSP', 'iEvalInfoMask');
        fMeasInfo(sFilename);
    elseif (lVer == 0)   % or dMDH(1,6) ?, bitget(uint16(dMDH(1, 1)), 16) == 0   %
        % DrecksMDH 1.0
        save(sFilenameMat, 'iLC', 'iSP', 'iEvalInfoMask', 'iLCDrecks', 'iLCPositions');
        fMeasInfo(sFilename);
    else
        % DrecksMDH 2.0
        dLCPositions = dMDH(~lDrecksMDH, 31:end);
        SPos.dSliceDataPosition = dLCPositions(:,1:3);
        if SDrecksMDH.Seq.Is3D
            SPos.dSliceDataPosition = unique(SPos.dSliceDataPosition, 'rows');
        end
        SPos.dTablePosition     = unique(-((dLCPositions(:,4) +5)/10));
        SPos.dRotMatrix         = quat2dcm(dLCPositions(1,5:8));
        SPos.dRotMatrix(:,2:3)  = -SPos.dRotMatrix(:,2:3); % due to negative orientation of y and z axis (compared to standard cartesian system)
        [yaw, pitch, roll]      = dcm2angle(SPos.dRotMatrix);
        SPos.dRotAngles         = [yaw, pitch, roll];  
        
        save(sFilenameMat, 'iLC', 'iSP', 'iEvalInfoMask', 'SDrecksMDH', 'SPos');
        fMeasInfo2(sFilename);       % Print loop counter summary
    end
else
    % No input argument: Execute on all *.dat files in folder (batch processing)
    fprintf(1, '\n***********************\n');
    fprintf(1, '* Entering batch-mode *\n');
    fprintf(1, '*   Get some coffee!  *\n');
    fprintf(1, '***********************\n\n');
    SDir = dir('*.dat');
    for iI = 1:length(SDir)
        fMeasCreateLUT_main(SDir(iI).name);
    end
end