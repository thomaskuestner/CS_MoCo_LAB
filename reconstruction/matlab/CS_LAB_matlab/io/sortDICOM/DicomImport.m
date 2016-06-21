function SDataOut = DicomImport(sPath, lRecursive, csAdditionalTags)
% DICOMIMPORT Import DICOM series with preview.
%
%   SDATA = DICOMIMPORT(sPATH, lRECURSIVE) Scans the folder specified by
%   sPATH for DICOM files. If lRECURSIVE is true, DICOMIMPORT
%   hierarchically scans all subolders of sPATH. Valid DICOM files are
%   sorted into series according the DICOM tag SeriesInstanceUID and the
%   folder in which the files were found. For each series, a thumbnail
%   image is created and an overview of the obtained series is displayed in
%   a GUI that lets you select the series to load. Select series as if you
%   would in a file browser by klicking and using the SHIFT and CNTL keys
%   to select multiple series Close the GUI by pressing <ENTER> or <ESCAPE>
%   (the latter returns an empty array). Output SDATA is a struct where
%   length(SDATA) equals the amout of selected series. SDATA contains at
%   least the following fields:
%
%       Img:                The image volume
%       SeriesDescription:  The contents of the corresponding DICOM tag
%       Orientation:        The most prominent image volume orientation
%       Aspect:             The pixel spacing in all 3 directions
%       ImageOrientation:   The direction cosines (compare DICOM standard)
%       ImagePosition:      The coordinates of the upper left corner of
%                           each image slice (compare DICOM standard)
%
%   The tags specified in csADDITIONALTAGS (see below) are appended to this
%   structure.
%
%
%   DICOMIMPORT lets you select the base folder in a dialog and asks if the
%   sobfolders are to be scanned.
%
%   DICOMIMPORT(sPATH, lRECURSIVE, csADDITIONALTAGS) Lets you specify
%   additional DICOM tags you can use for further processing in the cell
%   array of strings csADDITIONALTAGS
%
%
% Copyright 2013 Christian Wuerslin, University of Tuebingen and
% University of Stuttgart, Germany.
% Contact: christian.wuerslin@med.uni-tuebingen.de

% =========================================================================
% *** FUNCTION DicomImport
% ***
% *** Main GUI function. See above for description.
% ***
% =========================================================================

% -------------------------------------------------------------------------
% Control the figure's appearence
SAp.sTITLE              = 'Dicom Import [press enter to confirm]';
SAp.iTHUMBNAILSIZE      = 200;
SAp.iINFOBARHEIGHT      = 125;
SAp.dBGCOLOR            = [0.1 0.2 0.3];
SAp.iNROWS              = 2;
SAp.iNCOLS              = 4;
SAp.iTITLEHEIGHT        = 16;
SAp.dOPACITY            = 0.5;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define a list of DICOM tags to read (those that I really need)
csTags = {'SeriesInstanceUID',       'Rows',                 'Columns', ...
          'ImageOrientationPatient', 'ImagePositionPatient', 'PixelSpacing', ...
          'SeriesDescription', 'Modality'};
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parse input arguments
if nargin
    if ~ischar(sPath), error('Input argument must be a string!'); end
else
    sPath = uigetdir;
    if isnumeric(sPath)
        SDataOut = [];
        return
    end
end

if nargin < 2, lRecursive = strcmp(questdlg('Search subdirectories for DICOM files?', 'Stupid question', 'Yup', 'Noop', 'Yup'), 'Yup'); end
if nargin > 2, csTags = [csTags, csAdditionalTags]; end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Get all DICOM files and their headers
fprintf(1, 'Scanning for DICOM files...');
SData = fGetDicomFiles(sPath, csTags, lRecursive);
fprintf(1, 'done.\n');
if isempty(SData), return, end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Sort the series into a struct and get a thumbnail of each set
fprintf(1, 'Creating thumbnails...');
iHashs = [SData.iHash];
iSeriesHash = fUnique(iHashs);
iNDatasets = length(iSeriesHash);
iThumbnails = zeros(SAp.iTHUMBNAILSIZE, SAp.iTHUMBNAILSIZE, iNDatasets);
for iDataset = 1:iNDatasets
    iThisSeries = find([SData.iHash] == iSeriesHash(iDataset));
    SSeries(iDataset).SData = SData(iThisSeries);

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Load the center image and create thumbnail
    iInd = iThisSeries(round(length(iThisSeries)./2));
    dImg = double(dicomread(SData(iInd).sFilename));
    if isempty(dImg), continue, end
    
    dImg = dImg - min(dImg(:));
    dX = linspace(1, length(dImg), SAp.iTHUMBNAILSIZE);
    [dYY, dXX] = meshgrid(dX, dX);
    dImg = interp2(dImg, dYY, dXX, 'linear*', 0);
    if max(dImg(:)), dImg = dImg./max(dImg(:)); end
    iThumbnails(:,:,iDataset) = dImg;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end
lSelected = false(iNDatasets, 1);
fprintf(1, 'done.\n');
fprintf(1, 'Found %u files in %u datasets!\n', length(SData), iNDatasets);
clear iDataset iThisSeries iInd dImg dX dXX dYY
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Create the figure.
dFigureWidth  = (SAp.iTHUMBNAILSIZE + 1).*SAp.iNCOLS;
dFigureHeight = (SAp.iTHUMBNAILSIZE + SAp.iTITLEHEIGHT + 1).*SAp.iNROWS + SAp.iINFOBARHEIGHT + 2;
dScreenSize = get(0, 'MonitorPositions');
try
    hF = figure(...
        'Position'             , [(dScreenSize(3) - dFigureWidth )./2, ...
        (dScreenSize(4) - dFigureHeight)./2, ...
        dFigureWidth, dFigureHeight], ...
        'Units'                , 'pixels', ...
        'Color'                , SAp.dBGCOLOR, ...
        'Resize'               , 'off', ...
        'DockControls'         , 'off', ...
        'MenuBar'              , 'none', ...
        'Name'                 , SAp.sTITLE, ...
        'NumberTitle'          , 'off', ...
        'BusyAction'           , 'cancel', ...
        'KeyPressFcn'          , @fKeyPressFcn, ...
        'CloseRequestFcn'      , @fCloseGUI, ...
        'WindowButtonDownFcn'  , @fWindowButtonDownFcn, ...
        'WindowButtonMotionFcn', @fWindowMouseMoveFcn, ...
        'WindowScrollWheelFcn' , @fWindowScrollWheelFcn);
catch %#ok<CTCH>
    hF = figure(...
        'Position'             , [(dScreenSize(3) - dFigureWidth )./2, ...
        (dScreenSize(4) - dFigureHeight)./2, ...
        dFigureWidth, dFigureHeight], ...
        'Units'                , 'pixels', ...
        'Color'                , SAp.dBGCOLOR, ...
        'Resize'               , 'off', ...
        'DockControls'         , 'off', ...
        'MenuBar'              , 'none', ...
        'Name'                 , SAp.sTITLE, ...
        'NumberTitle'          , 'off', ...
        'BusyAction'           , 'cancel', ...
        'KeyPressFcn'          , @fKeyPressFcn, ...
        'CloseRequestFcn'      , @fCloseGUI, ...
        'WindowButtonDownFcn'  , @fWindowButtonDownFcn, ...
        'WindowButtonMotionFcn', @fWindowMouseMoveFcn);
end
clear dScreenSize
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Create the axes
hAxes = zeros(SAp.iNCOLS.*SAp.iNROWS, 1);
hImg  = zeros(SAp.iNCOLS.*SAp.iNROWS, 1);
hText = zeros(SAp.iNCOLS.*SAp.iNROWS, 1);
for iM = 1:SAp.iNROWS
    for iN = 1:SAp.iNCOLS
        iLinInd = (iM - 1).*SAp.iNCOLS + iN;
        dVertPos = dFigureHeight - iM.*(SAp.iTHUMBNAILSIZE + SAp.iTITLEHEIGHT + 1) + 1;
        iPosition = [(iN - 1).*(SAp.iTHUMBNAILSIZE + 1) + 2, dVertPos, SAp.iTHUMBNAILSIZE, SAp.iTHUMBNAILSIZE];
        hAxes(iLinInd) = axes('Parent' , hF, 'Units', 'pixels', 'Position', iPosition);
        hImg(iLinInd) = image(zeros(SAp.iTHUMBNAILSIZE, SAp.iTHUMBNAILSIZE, 3), 'Parent', hAxes(iLinInd));
        hText(iLinInd) = uicontrol(...
            'Parent'            , hF, ...
            'Style'             , 'text', ...
            'Units'             , 'pixels', ...
            'Position'          , [iPosition(1), iPosition(2) + iPosition(4), SAp.iTHUMBNAILSIZE, SAp.iTITLEHEIGHT], ...
            'BackgroundColor'   , SAp.dBGCOLOR./2, ...
            'ForegroundColor'   , 'w', ...
            'HorizontalAlign'   , 'left', ...
            'FontUnits'         , 'pixels', ...
            'FontSize'          , 12);
    end
end
axis(hAxes, 'off');
clear iM iN
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Create the info text objects
uicontrol(...
    'Parent'            , hF, ...
    'Style'             , 'text', ...
    'Units'             , 'pixels', ...
    'Position'          , [2, 2, 200, SAp.iINFOBARHEIGHT], ...
    'BackgroundColor'   , SAp.dBGCOLOR./2, ...
    'ForegroundColor'   , 'w', ...
    'HorizontalAlign'   , 'left', ...
    'FontUnits'         , 'pixels', ...
    'FontSize'          , 14, ...
    'String'            , {'Series Description:', 'Path:', 'Modality:', 'Rows:', 'Columns:', 'Images:', 'Pixel Spacing:'});
hInfoText = uicontrol(...
    'Parent'            , hF, ...
    'Style'             , 'text', ...
    'Units'             , 'pixels', ...
    'Position'          , [202, 2, dFigureWidth - 202, SAp.iINFOBARHEIGHT], ...
    'BackgroundColor'   , SAp.dBGCOLOR./2, ...
    'ForegroundColor'   , 'w', ...
    'HorizontalAlign'   , 'left', ...
    'FontUnits'         , 'pixels', ...
    'FontSize'          , 14, ...
    'String'            , '');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Create a beatuiful mask for the selection highlight
dMask = repmat(permute(SAp.dBGCOLOR, [1 3 2]) , [SAp.iTHUMBNAILSIZE, SAp.iTHUMBNAILSIZE, 1]);
dMask = dMask.*repmat(linspace(1, 0.5, SAp.iTHUMBNAILSIZE)', [1, SAp.iTHUMBNAILSIZE, 3]);
dMask = dMask + 0.05.*rand(size(dMask));
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Fill the figure
iStartSeries = 1;
iLastSeries = 0;
sResult = 'esc';
fFillPanels();
% -------------------------------------------------------------------------

uiwait(hF);
% -------------------------------------------------------------------------
% ~~~~~~~~~~~~ A lot of user-friendly GUI interaction going on ~~~~~~~~~~~~
% -------------------------------------------------------------------------


% The GUI was closed in some way, continue execution


% -------------------------------------------------------------------------
% Dialog was aborted: Seriously, who does that?
if strcmp(sResult, 'esc')
    SDataOut = [];
    try delete(hF); end
    return
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% No seriers selected? - Return
iInd = find(lSelected);
if isempty(iInd)
    SDataOut = [];
    delete(hF);
    return
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Load the selected series
fprintf(1, 'Loading slected DICOM files...');
for i = 1:length(iInd)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Load the specified images
    SThisData = SSeries(iInd(i)).SData;
    iNImages = length(SThisData);
    iImg = zeros(SThisData(1).Rows, SThisData(1).Columns, iNImages, 'uint16');
    for j = 1:iNImages
        iImg(:,:,j) = dicomread(SThisData(j).sFilename);
    end
    % Now images are in iImg, headers in SThisData
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Find out the major orientation
    dImageOrientation = reshape(SThisData(1).ImageOrientationPatient, [3, 2])'; % Should be the same for a volume
    dImagePosition = [SThisData.ImagePositionPatient]';
    
    dOrientIndicator = sum(abs(dImageOrientation));
    [temp, iMinInd] = min(dOrientIndicator); %#ok<ASGLU>
    d3rdDimInd = dImagePosition(:, iMinInd);
    [temp, iSortInd] = sort(d3rdDimInd, 'ascend'); %#ok<ASGLU>
    iImg = iImg(:,:,iSortInd);
    iImageOrientation = round(dImageOrientation);
    switch(iMinInd)
        case 1 % The x coordinate changes least in both directions -> Sagittal
            SDataOut(i).Orientation = 'Sag'; % -> x: a->p=pos, y: h->f=neg
            if iImageOrientation(1, 2) == -1, flipdim(iImg, 2); end % <- y-coordinate was flipped
            if iImageOrientation(2, 3) ==  1, flipdim(iImg, 1); end % <- z-coordinate was flipped
            
        case 2 % The y coordinate changes least in both directions -> Coronal
            SDataOut(i).Orientation = 'Cor'; % -> x: r->l=pos, y: h->f=neg
            if iImageOrientation(1, 1) == -1, flipdim(iImg, 2); end % <- x-coordinate was flipped
            if iImageOrientation(2, 3) ==  1, flipdim(iImg, 1); end % <- z-coordinate was flipped
            
        case 3 % The z coordinate changes least in both directions -> Transversal
            SDataOut(i).Orientation = 'Tra'; % -> x: r->l=pos, y: a->p=pos
            if iImageOrientation(1, 1) == -1, flipdim(iImg, 2); end % <- x-coordinate was flipped
            if iImageOrientation(2, 2) == -1, flipdim(iImg, 1); end % <- Y-coordinate was flipped
    end
    SDataOut(i).Img = iImg;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Save the most important infos to the ouput variable
    SDataOut(i).SeriesDescriptions = SThisData(1).SeriesDescription;
    SDataOut(i).Aspect = zeros(3, 1);
    SDataOut(i).Aspect(1:2) = SThisData(1).PixelSpacing;
    if ~isscalar(d3rdDimInd), SDataOut(i).Aspect(3)   = abs(d3rdDimInd(2) - d3rdDimInd(1)); end
    SDataOut(i).ImageOrientation = dImageOrientation;
    SDataOut(i).ImagePosition = dImagePosition(iSortInd, :);
    if nargin == 3
        for j = 1:length(csAdditionalTags) % The additionally requested infos
            sTag = csAdditionalTags{j};
            eval(['xData = [SThisData.', sTag, ']'';']);
            xData = xData(iSortInd);
            eval(['SDataOut(i).', sTag, ' = xData;']);
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end
% -------------------------------------------------------------------------
fprintf(1, 'done.\n');

delete(hF);
% =========================================================================
% ***
% *** The 'end' of the DICOMIMPORT main function.
% ***
% =========================================================================


 
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * NESTED FUNCTION fCloseGUI (nested in DicomImport)
    % * * 
    % * * Figure callback
    % * *
    % * * Closes the figure.
    % * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function fCloseGUI(hObject, eventdata) %#ok<*INUSD> eventdata is repeatedly unused
        uiresume(hObject);
        delete(hObject);
    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fCloseGUI
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * NESTED FUNCTION fWindowMouseMoveFcn (nested in DicomImport)
    % * * 
    % * * Figure callback
    % * *
    % * * Displays informations about the series under the mouse cursor in
    % * * the texts at the bottom of the figure.
    % * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function fWindowMouseMoveFcn(hObject, eventdata)
        if ~exist('hImg', 'var'), return, end; % Return if called during GUI startup
        
        iAxisInd = fGetAxes();
        iSeriesInd = iAxisInd + iStartSeries -1;
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Return if not over axes or out of data bounds
        if ~iAxisInd, set(hInfoText, 'String', ''); return, end
        if iSeriesInd > length(lSelected),  set(hInfoText, 'String', ''); return, end
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Print Information
        sText = SSeries(iSeriesInd).SData(1).SeriesDescription;
        if isempty(sText), sText = 'N/A'; end, csText(1) = {sText};
        
        sText = fileparts(SSeries(iSeriesInd).SData(1).sFilename);
        csText(2) = {sText};
        
        sText = SSeries(iSeriesInd).SData(1).Modality;
        if isempty(sText), sText = 'N/A'; end, csText(3) = {sText};
        
        sText = SSeries(iSeriesInd).SData(1).Rows;
        if isempty(sText), sText = 'N/A'; end, csText(4) = {sText};
        
        sText = SSeries(iSeriesInd).SData(1).Columns;
        if isempty(sText), sText = 'N/A'; end, csText(5) = {sText};
        
        sText = num2str(length(SSeries(iSeriesInd).SData));
        if isempty(sText), sText = 'N/A'; end, csText(6) = {sText};
        
        dPixelSpacing =  SSeries(iSeriesInd).SData(1).PixelSpacing;
        if ~isempty(dPixelSpacing)
            sText = sprintf('%2.2f x %2.2f mm', dPixelSpacing(1), dPixelSpacing(2));
        else
            sText = 'N/A';
        end
        csText(7) = {sText};
        
        set(hInfoText, 'String', csText);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fWindowMouseMoveFcn
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * NESTED FUNCTION fWindowButtonDownFcn (nested in DicomImport)
    % * * 
    % * * Figure callback
    % * *
    % * * Starting callback for mouse button actions. Manage the selection
    % * * of the series.
    % * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function fWindowButtonDownFcn(hObject, eventdata)
        iPanelInd = fGetAxes();
        iSeriesInd = iPanelInd + iStartSeries - 1;
        
        if ~iPanelInd, return, end % Exit if Event didn't occurr in axes
        if iSeriesInd > length(lSelected), return, end

        switch get(hF, 'SelectionType')
            
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Left mouse button without any key
            case 'normal'
                iNSelected = sum(lSelected);
                lState = lSelected(iSeriesInd);
                lSelected = false(size(lSelected));
                lSelected(iSeriesInd) = ~lState || iNSelected > 1;
                if ~lState, iLastSeries = iSeriesInd; end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Right mouse button/shift key
            case 'extend'
                if iLastSeries
                    iInd = sort([iLastSeries, iSeriesInd], 'ascend');
                    lSelected(iInd(1):iInd(2)) = true;
                else
                    lSelected(iSeriesInd) = true;
                    iLastSeries = iSeriesInd;
                end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Middle mouse button/alt key
            case 'alt'
                lState = lSelected(iSeriesInd);
                lSelected(iSeriesInd) = ~lState;
                if ~lState, iLastSeries = iSeriesInd; end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                
        end
        fFillPanels;
    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fWindowButtonDownFcn
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * NESTED FUNCTION fKeyPressFcn (nested in DicomImport)
    % * * 
    % * * Figure callback
    % * *
    % * * Callback for keyboard actions.
    % * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function fKeyPressFcn(hObject, eventdata) %#ok<INUSL>
    
        switch eventdata.Key

            case 'uparrow'
                iStartSeries = max([1, iStartSeries - SAp.iNCOLS]);
                fFillPanels;

            case 'downarrow'
                if iStartSeries + SAp.iNCOLS < length(lSelected),
                    iStartSeries = iStartSeries + SAp.iNCOLS;
                end
                
            case 'return' % Return and load the selected series
                sResult = 'OK';
                uiresume(hF);
                    
            case 'escape' % Return and do nothing
                uiresume(hF);
                
        end % switch
        
        fFillPanels;
        fWindowMouseMoveFcn(hF, []); % Make sure the info at the bottom is valid

    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fKeyPressFcn
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
       
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * NESTED FUNCTION fWindowScrollWheelFcn (nested in DicomImport)
    % * * 
    % * * Figure callback
    % * *
    % * * Callback for keyboard actions.
    % * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function fWindowScrollWheelFcn(hObject, eventdata)
        if eventdata.VerticalScrollCount > 0
            SEData.Key = 'downarrow';
        else
            SEData.Key = 'uparrow';
        end
        fKeyPressFcn(hObject, SEData); % This is the laszy way: just call the other callback
    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fWindowScrollWheelFcn
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    
       
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * *
    % * * NESTED FUNCTION fFillPanels (nested in DicomImport)
    % * *
    % * * Display the current data in all axes.
	% * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function fFillPanels()
        for iI = 1:length(hImg);
            iSeriesInd = iI + iStartSeries - 1;
            if iSeriesInd > length(lSelected)
                set(hImg(iI), 'CData', dMask./3);
                set(hText(iI), 'String', '<no data>');
            else
                dImg = repmat(iThumbnails(:,:,iSeriesInd), [1 1 3]);
                if lSelected(iSeriesInd),  dImg = 1 - ((1 - dImg).*(1 - SAp.dOPACITY.*dMask)); end
                set(hImg(iI), 'CData', dImg);
                sTitle = SSeries(iSeriesInd).SData(1).SeriesDescription;
                if length(sTitle) > 22, sTitle = [sTitle(1:21), '...']; end
                if isempty(sTitle), sTitle = '<no description>'; end
                set(hText(iI), 'String', sprintf('[%u]: %s', iSeriesInd, sTitle));
            end
        end
        drawnow expose
    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fFillPanels
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * *
    % * * NESTED FUNCTION fGetAxes (nested in DicomImport)
    % * *
    % * * Determine the axes number under the mouse cursor. Returns 0 if
    % * * not over an axes at all.
	% * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function iPanelInd = fGetAxes()
        iCursorPos = get(hF, 'CurrentPoint');
        iPanelInd = uint8(0);
        for iI = 1:length(hAxes)
            dPos = get(hAxes(iI), 'Position');
            if ((iCursorPos(1) >= dPos(1)) && (iCursorPos(1) < dPos(1) + dPos(3)) && ...
                (iCursorPos(2) >= dPos(2)) && (iCursorPos(2) < dPos(2) + dPos(4)))
                iPanelInd = uint8(iI);
            end
        end
    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fGetAxes
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * *
    % * * NESTED FUNCTION fGetDicomFiles (nested in DicomImport)
    % * *
    % * * Finds and creates a database of DICOM files
	% * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function SD = fGetDicomFiles(sFolder, csTags, lRecursive)
        SD = [];
        SDir = dir(sFolder);
        for iI = 1:length(SDir)
            sName = SDir(iI).name;
                        
            if (SDir(iI).isdir) && (lRecursive)
                if sName(1) == '.', continue, end
                SD = [SD, fGetDicomFiles([sFolder, filesep, SDir(iI).name], csTags, lRecursive)];
            else
                sFilename = [sFolder, filesep, SDir(iI).name];
                try
                    SHeader = dicominfo(sFilename);
                catch %#ok<CTCH>
                    continue
                end
                SRecord = [];
                for iJ = 1:length(csTags)
                    sTag = csTags{iJ};
                    try
                        eval(['SRecord.', sTag, ' = SHeader.', sTag, ';']);
                    catch %#ok<CTCH>
                        eval(['SRecord.', sTag, ' = [];']);
                    end
                end
                if isfield(SHeader, 'SeriesInstanceUID')
                    SRecord.iHash = fHash([sFolder, SHeader.SeriesInstanceUID]);
                    SRecord.sFilename = sFilename;
                    SD = [SD, SRecord];
                end
            end
        end
    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fGetDicomFiles
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * *
    % * * NESTED FUNCTION fHash (nested in DicomImport)
    % * *
    % * * Return a hash number
	% * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function iHashVal = fHash(sString)
        dString = double(sString);
        dHashVal = 5381;
        for iI = 1:length(dString)
            dHashVal = mod(dHashVal * 33 + dString(iI), 2.^32 - 1);
        end
        iHashVal = uint32(dHashVal);
    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fHash
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * *
    % * * NESTED FUNCTION fUnique (nested in DicomImport)
    % * *
    % * * Returns unique values in order of appearance
	% * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function xUniqueVals = fUnique(xArray)
        xUniqueVals = zeros(length(xArray), 1, class(xArray));
        iInd = 0;
        while ~isempty(xArray)
            iInd = iInd + 1;
            xUniqueVals(iInd) = xArray(1);
            xArray(xArray == xArray(1)) = [];
        end
        xUniqueVals = xUniqueVals(1:iInd);
    end
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fUnique
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

end
