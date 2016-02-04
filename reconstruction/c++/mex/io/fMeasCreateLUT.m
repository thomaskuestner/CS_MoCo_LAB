function fMeasCreateLUT(sFilename)
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
%   Copyright 2012 Christian Wuerslin, University of Tuebingen, Germany
%   christian.wuerslin@med.uni-tuebingen.de
%   Modifications 2013 Thomas Kuestner, University of Tuebingen, Germany
%   thomas.kuestner@med.uni-tuebingen.de

iMAXSIZE = 6000000;                     % ~200 MB memory usage
% This is the maximum allowed number of ADCs in your meas-file. Memory
% usage is 40*iMAXSIZE bytes. Adjust to your needs.

% csLCS = {'Smp', 'Cha', 'Lin', 'Acq', 'Sli', 'Par', 'Eco', 'Pha', 'Rep', 'Set', 'Seg', 'FreePara1', 'FreePara2', 'FreePara3', 'FreePara4'};
% Loop Counter names

if nargin
    % parse specified file
    SDir = dir(sFilename);
    if isempty(SDir)
        error('File doesn''t exist!');
    end
    iFileSize = SDir(1).bytes;
    fid = fopen(sFilename, 'r');
    iHeaderLength = fread(fid, 1, 'uint32');            % Read Header Length
    fprintf(1, '\n*** Parsing %s ***\n', sFilename);
    fprintf(1, 'File Size            : %u  (%u MB)\n', iFileSize, round(iFileSize/(1024^2)));
    fprintf(1, 'Length of the Header : %u\n', iHeaderLength);
    fseek(fid, iHeaderLength, 'bof');                   % Skip to first MDH
    iPos = ftell(fid);

    lLastADC = false;
    iCounter = 0;
    iCounterDrecks = 0;
    iLC = zeros(iMAXSIZE, 15, 'uint16');                 % Will hold No. of samples, No. of channels and loop counters
    iLCPositions = zeros(iMAXSIZE, 8, 'single');         % Will hold the slice data position vector (sagittal, coronal, transversal shift [mm] out of isocenter and table position [mm])
    iSP = zeros(iMAXSIZE,  1, 'uint64');                 % Will hold start of each MDH
    iLCDrecks = zeros(20, 41, 'uint16');                 % Will hold DrecksMDH values
    iEvalInfoMask = zeros(iMAXSIZE, 2,'uint32');

    fprintf(1, 'Parsing measurement file...');

    % Start the "Scout": Get No of MDHs, No of samples, No of channels and loop counters for each MDH
    while ~feof(fid) && ~lLastADC
        iCounter = iCounter + 1;
        iSP(iCounter) = iPos;                              % Start index of MDH
        fseek(fid, 20, 0);
        iEvalInfoMask(iCounter, :) = fread(fid, 2, 'uint32');
%         fseek(fid, 28, 0);                                 % skip to no of samples in scan
%         iLC(iCounter, :) = fread(fid, 16, 'uint16');
%         tmpLC = fread(fid, 40, 'uint16');
        tmpLC = [fread(fid, 20, 'uint16'); fread(fid, 1, 'float'); fread(fid, 1, 'uint32'); fread(fid, 10, 'uint16'); fread(fid, 7, 'float'); fread(fid, 2, 'uint16')];
        if(tmpLC(1) ~= 20)
            iLC(iCounter, :) = [tmpLC(1:11).', tmpLC(29:32).'];       % Read No of samples, No of channels and loop counters / for consistency -> leave it 16 samples long
            iLCPositions(iCounter, :) = [tmpLC(33:35).', tmpLC(41), tmpLC(36:39).'];
            fseek(fid, double(iLC(iCounter, 1)*8), 0);      % total MDH length = 128 Bytes
        else
            iCounterDrecks = iCounterDrecks + 1;
            iLCDrecks(iCounterDrecks,:) = tmpLC.';
            fseek(fid, double(iLCDrecks(iCounterDrecks,1)*8), 0);
        end
%         fseek(fid, double(68 + iLC(iCounter, 1)*8), 0);            % skip to beginning of next MDH (68 Bytes are the remainder of the MDH, iLC(...)*8 is the length of the ADC)
%         fseek(fid, double(20 + iLC(iCounter, 1)*8), 0);
        
        iPos = ftell(fid);                                 % Get index
        if iPos >= iFileSize - 384, lLastADC = true; end   % exclude the last 384 bytes of the file
        if ~mod(iCounter, 10000), fprintf(1, '.'); end     % show that something is going on
        if(iCounter > iMAXSIZE), error('Measurement file too large, increase iMAXSIZE'); end
    end

    fclose(fid);
    fprintf(1, 'done.\n');
    fprintf(1, 'Number of ADCs        : %u\n', iCounter);

    iLC = iLC(1:iCounter, :);   % Crop loop counters matrix
    iLCPositions = iLCPositions(1:iCounter,:);   % Crop position loop counters matrix
    iSP = iSP(1:iCounter);      % Crop start positions array
    iLCDrecks = iLCDrecks(1:iCounterDrecks,:); % Crop dreck loop counters matrix
    iEvalInfoMask = iEvalInfoMask(1:iCounter,:); % Crop Eval info mask

    [sPath, sName] = fileparts(sFilename);
    if isempty(sPath)
        save([sName, '.mat'], 'iLC', 'iSP', 'iLCDrecks', 'iLCPositions', 'iEvalInfoMask');
    else
        save([sPath, filesep, sName, '.mat'], 'iLC', 'iSP', 'iLCDrecks', 'iLCPositions', 'iEvalInfoMask');
    end
    
    fMeasInfo(sFilename);       % Print loop counter summary
    
else
    % No input argument: Execute on all *.dat files in folder (batch processing)
    fprintf(1, '\n***********************\n');
    fprintf(1, '* Entering batch-mode *\n');
    fprintf(1, '*   Get some coffee!  *\n');
    fprintf(1, '***********************\n\n');
    SDir = dir('*.dat');
    for iI = 1:length(SDir)
        fMeasCreateLUT(SDir(iI).name);
    end
end