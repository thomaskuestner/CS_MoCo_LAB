function [dData, iLC, iEvalInfoMask] = fMeasRead(sFilename, varargin)
%FMEASREAD Reads ADC data from Siemens measdata file.
%   [DDATA, ILC] = FMEASREAD(SFILENAME) tries to load all MDHs from the
%   Siemens measdata file SFILENAME provided the corresponding look-up
%   table file was created using FMEASCREATELUT. Output DDATA contains the
%   complex ADC data and ILC the loop counters of the actually loaded MDHs
%   for sorting.
% 
%   NOTE: fMeasRead returns an error if the length of the selected ADCs
%   differs. For some reason, the scanner usually adds an MDH at the end of
%   the measurement with length 16 or something. Thus calling fMeasRead
%   without a restriction on 'Smp' will usually cause an error. Try calling
%   fMeasRead(sFilename, 'Smp', [17 4096]) to get all measured ADCs.
%   
%   [ ... ] = FMEASREAD(SFILENAMES, 'PARAM1', val1, 'PARAM2', val2', ...)
%   specifies optional parameter name/value pairs to restrict the range of
%   an arbitrary selection of loop counters to a specific value/range.
%   Values can be either scalars, restricting the loop counter to one value
%   or two-element vectors specifying a range. Parameters/loop counters are:
%
%   'Smp' - No. of samples in ADC
%   'Cha' - No. of used channels
%   'Lin' - line LC
%   'Acq' - Acquisition LC
%   'Sli' - Slice LC
%   'Par' - Partition LC
%   'Eco' - Echo LC
%   'Pha' - Phase LC
%   'Rep' - Repetition LC
%   'Set' - Set LC
%   'Seg' - Segment LC
%   'FreePara1' - Free parameter 1 LC
%   'FreePara2' - Free parameter 2 LC
%   'FreePara3' - Free parameter 3 LC
%   'FreePara4' - Free parameter 4 LC
%
%   Example:
%       [dData, iLC] = fMeasRead('meas.dat', 'Sli', 1, 'Rep', [1 10]);
%       Reads repetions 1 to 10 from slice 1 in meas.dat.
%
%   Class support for input sFilename:
%      string
%
%   Class support for input varargin:
%       Pairs of string and integer
%
%   NOTE: All fMeas... routines assume that the channels is the innermost
%   dimension in the measdata file (which is true for all files I've seen
%   so far. Thus, when sorting the data into separate channels for
%   reconstruction, use the colon operator in the fashion
%   chan1data = alldata(1:NoChans:end, :);
%   chan2data = alldata(2:NoChans:end, :);
%   and so on.
%
%   See also FMEASCREATELUT, FMEASINFO
%
%   Copyright 2012 Christian Wuerslin, University of Tuebingen, Germany
%   christian.wuerslin@med.uni-tuebingen.de
%   Modifications 2013 Thomas Kuestner, University of Tuebingen, Germany
%   thomas.kuestner@med.uni-tuebingen.de

csLCS = {'Smp', 'Cha', 'Lin', 'Acq', 'Sli', 'Par', 'Eco', 'Pha', 'Rep', 'Set', 'Seg', 'FreePara1', 'FreePara2', 'FreePara3', 'FreePara4'};

if nargin == 0, error('No filename given.'); end

% Load the loop counters
[sPath, sName] = fileparts(sFilename);
if isempty(sPath)
    load([sName, '.mat']);
else
    load([sPath, filesep, sName, '.mat']);
end

% Set default values
iMin = min(iLC);
iMax = max(iLC);

iNOptInp = size(varargin, 2);
if mod(iNOptInp, 2)
    error('"varargin" must contain pairs of strings and integers!');
end

% Incorporate changes due to optional inputs (loop counter restrictions)
for iParam = 1:2:iNOptInp
    sCounter = varargin{iParam};
    iLimits = varargin{iParam + 1};
    
    if ~ischar(sCounter), error('Varargin must be pairs of strings and values!'); end
    if ~isnumeric(iLimits), error('Varargin must be pairs of strings and values!'); end
    
    for iCounter = 1:length(csLCS)
        if ~strcmp(sCounter, csLCS{iCounter}), continue; end
        iMin(iCounter) = iLimits(1);
        iMax(iCounter) = iLimits(end);
    end
end

% Create mask
lMask = true(size(iLC, 1), 1);
for iCounter = 1:size(iMin, 2)
    lMask = lMask & ((iLC(:, iCounter) >= iMin(iCounter)) & (iLC(:, iCounter) <= iMax(iCounter)));
end

% Crop output LC and SP and check if data is consistent
iLC = iLC(lMask, :);
iSP = iSP(lMask);
if(exist('iEvalInfoMask','var'))
    iEvalInfoMask = iEvalInfoMask(lMask,:);
else
    iEvalInfoMask = [];
end
if max(iLC(:,1)) ~= min(iLC(:,1)) % possible for Navis
    warning('Number of sampling points in selected ADCs not constant!');
    % dirty workaround => no. samples are most occuring ones
    nSmp = unique(iLC(:,1));
    helper = zeros(1,length(nSmp));
    for i=1:length(nSmp)
        helper(i) = sum(iLC(:,1) == nSmp(i));
    end
    [~, helper] = max(helper);
    lMask = iLC(:,1) == nSmp(helper);
    iLC = iLC(lMask, :);
    iSP = iSP(lMask);
    if(~isempty(iEvalInfoMask))
        iEvalInfoMask = iEvalInfoMask(lMask, :);
    end
end
iNSamples = iLC(1, 1);


% Load the data
fprintf(1, 'Number of matching ADCs: %u\n', length(iSP));
dRealData = zeros(length(iSP), iNSamples,'single');
dImagData = zeros(length(iSP), iNSamples,'single');

if(usejava('jvm') && ~feature('ShowFigureWindows'))
    flagDisp = false;
else
    flagDisp = true;
    if(exist('multiWaitbar','file'))
        flagMW = true;
    else
        flagMW = false;
    end
end

if(flagDisp), if(flagMW), multiWaitbar('Loading specified ADCs', 0); else hw = waitbar(0,'Loading specified ADCs'); end; end;
fid = fopen([sPath, filesep, sName, '.dat'], 'r');
for iI = 1:length(iSP)
    fseek(fid, double(iSP(iI)) + 128, 'bof');
    dLine = fread(fid, double(iNSamples*2), 'float');
    dRealData(iI, :) = dLine(1:2:end);
    dImagData(iI, :) = dLine(2:2:end);
    if(flagDisp), if(flagMW), multiWaitbar('Loading specified ADCs', iI/length(iSP)); else waitbar(iI/length(iSP),hw); end; end;
end
fclose(fid);
if(flagDisp), if(flagMW), multiWaitbar('Loading specified ADCs', 'Close'); else close(hw); end; end;

dData = complex(dRealData, dImagData);

















