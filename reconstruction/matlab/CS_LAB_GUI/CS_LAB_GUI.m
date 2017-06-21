function varargout = CS_LAB_GUI(varargin)
% CS_LAB_GUI MATLAB code for CS_LAB_GUI.fig
%      CS_LAB_GUI, by itself, creates a new CS_LAB_GUI or raises the existing
%      singleton*.
%
%      H = CS_LAB_GUI returns the handle to a new CS_LAB_GUI or the handle to
%      the existing singleton*.
%
%      CS_LAB_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CS_LAB_GUI.M with the given input arguments.
%
%      CS_LAB_GUI('Property','Value',...) creates a new CS_LAB_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CS_LAB_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CS_LAB_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CS_LAB_GUI

% Last Modified by GUIDE v2.5 15-Feb-2016 11:19:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CS_LAB_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CS_LAB_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CS_LAB_GUI is made visible.
function CS_LAB_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CS_LAB_GUI (see VARARGIN)

currpath = fileparts(mfilename('fullpath'));
cd(currpath);
addpath(genpath(currpath));
addpath(genpath([currpath,filesep,'utils']));
addpath(genpath([fileparts(currpath),filesep,'CS_LAB']));
addpath(genpath([fileparts(currpath),filesep,'CS_LAB_matlab']));
handles.currpath = currpath;

warning('off','MATLAB:mat2cell:TrailingUnityVectorArgRemoved');
warning('off','MATLAB:DELETE:FileNotFound');

% set default values
handles.sPathPDexe = 'sampling';
handles.cCSAlgos = {'FOCUSS', 'BART', 'sparseMRI', 'SPIRiT_CG', 'SPIRiT_POCS', 'ESPIRiT_CG', 'ESPIRiT_L1', 'L1_Magic_TV', 'L1_Magic_L1', 'L1_Magic_TVDantzig', 'L1_Magic_L1Dantzig', 'FCSA_WaTMRI', 'FCSA_SLEP', 'FCSA_proxA', 'BFCSA_proxA', 'ADMM_proxA', 'SB_proxA', 'Zero'};
handles.cSparseTrafos = {'fft', 'pca', 'dct', 'wavelet_mat', 'wavelet_lab', 'mellin', 'surfacelet'};
sFiles = dir([currpath,filesep,'datasets']);
handles.cDatasets = cell(0,0);
for iI = 1:length(sFiles)
    if(sFiles(iI).isdir), continue; end;
    [~, sFilename] = fileparts(sFiles(iI).name);
    handles.cDatasets(:,end+1) = {sFilename; [currpath,filesep,'datasets',filesep,sFiles(iI).name]};
end
% handles.cDatasets = {'brain 8ch (retrospective)', 'brain 1ch (retrospective)', 'phantom (retrospective)', 'resolution phantom (retrospective)', 'knee (prospective)', 'angio (prospective)', 'thorax 3x (prospective)', 'thorax 8x (prospective)', 'abdomen 4x (prospective)', 'abdomen 6x (prospective)', 'resolution phantom 4x (prospective)', 'resolution phantom 14x (prospective)'; ...
%                     'brain_full.mat', 'brain_1ch_full.mat', 'phantom_full.mat', 'phantom_1x_sub.mat', 'knee_sub.mat', 'angio_sub.mat', 'thorax_3x_sub.mat', 'thorax_8x_sub.mat', 'abdomen_4x_sub.mat', 'abdomen_6x_sub.mat', 'phantom_4x_sub.mat', 'phantom_14x_sub.mat'};
% all input data must be in a cell (1 x #Channels) with ky - kx - kz|slice
% => for CS_reconstruction assumed sparse dimensions are always 
handles.dRange = [0 1];
handles.cPossible = {true(1,length(handles.cSparseTrafos)), true(1,length(handles.cSparseTrafos))};
handles.idnames = {'axInput', 'axMask', 'axOut_1', 'axOut_2', 'axDiff_1', 'axDiff_2'; ...
    'dInput', 'dMask', 'dOut_1', 'dOut_2', 'dDiff_1', 'dDiff_2'};
handles.dimMask = [];
handles.iEsp = [1, 0];

% set popup menus
set(handles.popDatasets,'String', handles.cDatasets(1,:));
set(handles.popCSAlgo_1,'String', handles.cCSAlgos);
set(handles.popSparseTrafo_1, 'String', handles.cSparseTrafos);
set(handles.popCSAlgo_2,'String', handles.cCSAlgos);
set(handles.popSparseTrafo_2, 'String', handles.cSparseTrafos);

[handles.dKSpace, handles.dKSpaceShow, handles.dImg, handles.dMask, handles.measPara, sInfoTxt] = fLoadDataset(handles.cDatasets{2,1}, get(handles.edSparse,'Value'));
handles.iCurrSlice = round(size(handles.dImg,3)/2);
handles.dRecon_1 = [];
handles.dRecon_2 = [];
handles.dInput = scaleImg(handles.dImg, handles.dRange);
% initialize shown axes
defaultImg = zeros(200,200);
for iI=1:size(handles.idnames,2)
    if(iI > 2)
        handles.(handles.idnames{2,iI}) = defaultImg;
    end
    handles.hPlot.(handles.idnames{2,iI}) = fPlotAxes(handles.(handles.idnames{2,iI}), 1, [], handles.(handles.idnames{1,iI}), handles.dRange);
end
linkaxes([handles.axInput, handles.axOut_1, handles.axOut_2, handles.axDiff_1, handles.axDiff_2]);

set(handles.txtInfo,'String',sInfoTxt{1});
set(handles.txtMask,'String',sInfoTxt{2});
if(~isempty(handles.dMask)) % mask was supplied
    handles.measPara.prospectiveSub = true;
    handles.dimMask = size(handles.dMask);
    set(handles.edSparse,'Value',handles.measPara.iIndSparseDim);   
    set(handles.edSparse,'Enable','off');
else
    handles.measPara.prospectiveSub = false;
    set(handles.edSparse,'Enable','on');
end   
if(size(handles.dImg,3) == 1) % just 2D image
    set(handles.popTrafoDim_1,'Value',2);
    set(handles.popTrafoDim_2,'Value',2);
else % 3D image
    set(handles.popTrafoDim_1,'Value',3);
    set(handles.popTrafoDim_2,'Value',3);
end

% set icons
dImage = double(imread(['icons',filesep,'open.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pb_LoadData, 'CData', dImage/max(dImage(:)));
set(handles.pbLoad_1, 'CData', dImage/max(dImage(:)));
set(handles.pbLoad_2, 'CData', dImage/max(dImage(:)));

dImage = double(imread(['icons',filesep,'save.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pbSave_1, 'CData', dImage/max(dImage(:)));
set(handles.pbSave_2, 'CData', dImage/max(dImage(:)));

dImage = double(imread(['icons',filesep,'reset.png']))./255;
if size(dImage, 3) == 1, dImage = repmat(dImage, [1 1 3]); end
set(handles.pbReset, 'CData', dImage/max(dImage(:)));

% delete previously existing parameter files
delete('parameter_currRun_1.m', 'parameter_currRun_2.m', 'samplingPattern.txt');

% Choose default command line output for CS_LAB_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CS_LAB_GUI wait for user response (see UIRESUME)
% uiwait(handles.CS_LAB_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = CS_LAB_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function popDatasets_Callback(hObject, eventdata, handles)
% change datasets

if(nnz(handles.dRecon_1) > 0 || nnz(handles.dRecon_2) > 0)
    choice = questdlg('Existing reconstruction will be deleted. Do you want to continue?', 'Deleting results', 'Yes', 'No', 'Yes');
    switch choice
        case 'Yes'
            if(exist([handles.currpath,filesep,'parameter_currRun_1.m'],'file'))
                delete([handles.currpath,filesep,'parameter_currRun_1.m']);
            end
            if(exist([handles.currpath,filesep,'parameter_currRun_2.m'],'file'))
                delete([handles.currpath,filesep,'parameter_currRun_2.m']);
            end
        case 'No'
            return;
    end
end
[handles.dKSpace, handles.dKSpaceShow, handles.dImg, dMask, handles.measPara, sInfoTxt] = fLoadDataset(handles.cDatasets{2,get(hObject,'Value')}, get(handles.edSparse,'Value'));
if(~isempty(dMask))
    handles.measPara.prospectiveSub = true;    
    if(~isempty(handles.dMask) && all(size(dMask) == size(handles.dMask)))
        choice = questdlg('Subsampling mask already existing. Overwrite it with one from file?', 'Subsampling mask', 'Yes', 'No', 'No');
        switch choice
            case 'Yes'
                handles.dMask = dMask;
                handles.dimMask = size(handles.dMask);
                set(handles.edSparse,'Value',handles.measPara.iIndSparseDim);   
                set(handles.edSparse,'Enable','off');
            case 'No'
                set(handles.edSparse,'Enable','on');
                % NOP
        end
    else
        handles.dMask = dMask;
        handles.dimMask = size(handles.dMask);
        set(handles.edSparse,'Value',handles.measPara.iIndSparseDim);   
        set(handles.edSparse,'Enable','off');
    end
else
    handles.dMask = [];
    handles.measPara.prospectiveSub = false;
    set(handles.edSparse,'Enable','on');
end
if(size(handles.dImg,3) == 1) % just 2D image
    set(handles.popTrafoDim_1,'Value',2);
    set(handles.popTrafoDim_2,'Value',2);
else % 3D image
    set(handles.popTrafoDim_1,'Value',3);
    set(handles.popTrafoDim_2,'Value',3);
end
handles.dOut_1 = zeros(size(handles.dImg));
handles.dOut_2 = zeros(size(handles.dImg));
handles.dRecon_1 = [];
handles.dRecon_2 = [];
handles.dDiff_1 = zeros(size(handles.dImg));
handles.dDiff_2 = zeros(size(handles.dImg));

% get shown parameter
iInd = get(handles.popShowInput, 'Value');
switch iInd
    case 1 % image
        handles.dInput = scaleImg(handles.dImg, handles.dRange);        
    case 2 % kspace
        handles.dInput = scaleImg(handles.dKSpaceShow, handles.dRange);    
end
handles.iCurrSlice = round(size(handles.dInput,3)/2);
for iI=1:size(handles.idnames,2)
    handles.hPlot.(handles.idnames{2,iI}) = fPlotAxes(handles.(handles.idnames{2,iI}), handles.iCurrSlice, handles.hPlot.(handles.idnames{2,iI}), handles.(handles.idnames{1,iI}), handles.dRange);
end
set(handles.txtInfo,'String',sInfoTxt{1});
set(handles.txtMask,'String',sInfoTxt{2});
guidata(hObject, handles);        


% --- Executes during object creation, after setting all properties.
function popDatasets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popDatasets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popShowInput_Callback(hObject, eventdata, handles)
% change between image and kspace

iInd = get(hObject,'Value');

switch iInd
    case 1 % image
        handles.dInput = scaleImg(handles.dImg, handles.dRange);        
    case 2 % kspace
        handles.dInput = scaleImg(handles.dKSpaceShow, handles.dRange);    
end
handles.hPlot.dInput = fPlotAxes(handles.dInput, handles.iCurrSlice, handles.hPlot.dInput, [], []);
guidata(hObject, handles);
        


% --- Executes during object creation, after setting all properties.
function popShowInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popShowInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pbCreateMask_Callback(hObject, eventdata, handles)
% create subsampling mask

if(~isempty(handles.dMask))
    choice = questdlg('Subsampling mask already existing. Overwrite with new one?', 'Subsampling mask', 'Yes', 'No', 'Yes');
    switch choice
        case 'Yes'
            % NOP
        case 'No'
            return;
    end
end
cd(handles.currpath);
handles = fCreateMask(handles);
guidata(hObject,handles);

function handles = fCreateMask(handles)
iIndSparseDim = get(handles.edSparse,'Value');
switch iIndSparseDim
    case 1 % yPhase, zPhase
        dimMask = [handles.measPara.dim(1), handles.measPara.dim(3)];
    case 2 % yPhase, xFreq
        dimMask = [handles.measPara.dim(1), handles.measPara.dim(2)];
    case 3 % xFreq, zPhase
        dimMask = [handles.measPara.dim(2), handles.measPara.dim(3)];
    case 4 % yPhase
        dimMask = [handles.measPara.dim(1), 1];
    case 5 % xFreq
        dimMask = [handles.measPara.dim(2), 1];
    case 6 % zPhase
        dimMask = [handles.measPara.dim(3), 1];        
end
handles.dimMask = dimMask;
dAccel = str2double(get(handles.edAccel,'String'));
if(dAccel < 2 || dAccel > 15)
    vname = @(x) inputname(1);
    errordlg(sprintf('%s: value not accepted', vname(dAccel)));
end
dFullySampled = str2double(get(handles.edFullySampled,'String'));
if(dFullySampled < 0 || dFullySampled > 1)
    vname = @(x) inputname(1);
    errordlg(sprintf('%s: value not accepted', vname(dFullySampled)));
end
if(dFullySampled * prod(dimMask) > prod(dimMask)/dAccel)
    errordlg('number of fully_sampled-region-points is higher then the number of points to create!');
end
iVD = get(handles.popVD,'Value');
iSubAlgo = get(handles.popSubAlgo,'Value');
iEsp = get(handles.popESPReSSo,'Value');
switch iEsp
    case 1
        iEspIn = [1, 0]; % factor, direction
    case 2
        iEspIn = [7/8, 1];
    case 3
        iEspIn = [6/8, 1];
    case 4
        iEspIn = [5/8, 1];
    case 5
        iEspIn = [4/8, 1];
    case 6
        iEspIn = [7/8, 0];
    case 7
        iEspIn = [6/8, 0];
    case 8
        iEspIn = [5/8, 0];
    case 9
        iEspIn = [4/8, 0];
end
handles.iEsp = iEspIn;
iEllScanning = get(handles.cbEllScanning,'Value') + 1;
sTmp = {'false','true'};
if(ispc), sExt = '.exe'; else sExt = ''; end
system([handles.sPathPDexe,filesep,'Subsample',sExt,sprintf(' %d %d %f %d %f %s %.3f %s %d', dimMask(1), dimMask(2), dAccel, iSubAlgo, dFullySampled, sTmp{iEllScanning}, iEspIn(1), sTmp{iEspIn(2)+1}, iVD)]);
handles.dMask = readSamplingMaskFromFile('samplingPattern.txt');

if(size(handles.dMask,2) == 1)
    handles.dMask = repmat(handles.dMask,1,size(handles.dInput,2));
end
if(iIndSparseDim == 5) % xFreq
    handles.dMask = handles.dMask.';
end
handles.hPlot.dMask = fPlotAxes(handles.dMask, 1, handles.hPlot.dMask, handles.axMask, handles.dRange);


function popCSAlgo_1_Callback(hObject, eventdata, handles)
% change reconstruction algorithm of algo1

iInd = get(hObject,'Value');
switch handles.cCSAlgos{iInd}
    case {'FOCUSS', 'sparseMRI'}
        handles.cPossible{1} = true(1,length(handles.cSparseTrafos));
    case {'SPIRiT_CG', 'SPIRiT_POCS', 'ESPIRiT_CG', 'ESPIRiT_L1'}
        handles.cPossible{1} = logical([0 0 0 0 1 0 0]);
    case {'L1_Magic_TV', 'L1_Magic_L1', 'L1_Magic_TVDantzig', 'L1_Magic_L1Dantzig'}
        handles.cPossible{1} = logical([1 0 0 0 0 0 0]);
    case {'BART', 'FCSA_WaTMRI', 'FCSA_SLEP', 'FCSA_proxA', 'BFCSA_proxA', 'ADMM_proxA', 'SB_proxA'}
        handles.cPossible{1} = logical([0 0 0 1 0 0 0]);
    case 'Zero'
        handles.cPossible{1} = logical([1 0 0 0 0 0 0]);
end
if(exist([handles.currpath,filesep,'parameter_currRun_1.m'],'file'))
    delete([handles.currpath,filesep,'parameter_currRun_1.m']);
end
set(handles.popSparseTrafo_1, 'String', handles.cSparseTrafos(handles.cPossible{1}), 'Value', 1);
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function popCSAlgo_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popCSAlgo_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popSparseTrafo_1_Callback(hObject, eventdata, handles)
% change sparsifying transformation of algo1

% Hints: contents = cellstr(get(hObject,'String')) returns popSparseTrafo_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popSparseTrafo_1
if(exist([handles.currpath,filesep,'parameter_currRun_1.m'],'file'))
    delete([handles.currpath,filesep,'parameter_currRun_1.m']);
end

% --- Executes during object creation, after setting all properties.
function popSparseTrafo_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popSparseTrafo_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pbCSEdit_1_Callback(hObject, eventdata, handles)
% edit reconstruction parameter of algo1
if(~exist([handles.currpath,filesep,'parameter_currRun_1.m'],'file'))
    copyfile([handles.currpath,filesep,'parameters_default_GUI.m'], [handles.currpath,filesep,'parameter_currRun_1.m']);
end
contents = cellstr(get(handles.popCSAlgo_1,'String'));
sAlgo = contents{get(handles.popCSAlgo_1,'Value')};
fReplaceText([handles.currpath,filesep,'parameter_currRun_1.m'], sAlgo, 'cstype');
contents = cellstr(get(handles.popSparseTrafo_1,'String'));
sTrafo = contents{get(handles.popSparseTrafo_1,'Value')};
fReplaceText([handles.currpath,filesep,'parameter_currRun_1.m'], sTrafo, 'transformation');
edit([handles.currpath,filesep,'parameter_currRun_1.m']);


function popCSAlgo_2_Callback(hObject, eventdata, handles)
% change reconstruction algorithm of algo2

iInd = get(hObject,'Value');
switch handles.cCSAlgos{iInd}
    case {'FOCUSS', 'sparseMRI'}
        handles.cPossible{2} = true(1,length(handles.cSparseTrafos));
    case {'SPIRiT_CG', 'SPIRiT_POCS', 'ESPIRiT_CG', 'ESPIRiT_L1'}
        handles.cPossible{2} = logical([0 0 0 0 1 0 0]);
    case {'L1_Magic_TV', 'L1_Magic_L1', 'L1_Magic_TVDantzig', 'L1_Magic_L1Dantzig'}
        handles.cPossible{2} = logical([1 0 0 0 0 0 0]);
    case {'BART', 'FCSA_WaTMRI', 'FCSA_SLEP', 'FCSA_proxA', 'BFCSA_proxA', 'ADMM_proxA', 'SB_proxA'}
        handles.cPossible{2} = logical([0 0 0 1 0 0 0]);
    case 'Zero'
        handles.cPossible{2} = logical([1 0 0 0 0 0 0]);
end
if(exist([handles.currpath,filesep,'parameter_currRun_2.m'],'file'))
    delete([handles.currpath,filesep,'parameter_currRun_2.m']);
end
set(handles.popSparseTrafo_2, 'String', handles.cSparseTrafos(handles.cPossible{2}) ,'Value', 1);
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popCSAlgo_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popCSAlgo_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popSparseTrafo_2_Callback(hObject, eventdata, handles)
% change sparsifying transformation of algo2

if(exist([handles.currpath,filesep,'parameter_currRun_2.m'],'file'))
    delete([handles.currpath,filesep,'parameter_currRun_2.m']);
end


function popSparseTrafo_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popSparseTrafo_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pbCSEdit_2_Callback(hObject, eventdata, handles)
% edit reconstruction parameters of algo2
if(~exist([handles.currpath,filesep,'parameter_currRun_2.m'],'file'))
    copyfile([handles.currpath,filesep,'parameters_default_GUI.m'], [handles.currpath,filesep,'parameter_currRun_2.m']);
end
contents = cellstr(get(handles.popCSAlgo_2,'String'));
sAlgo = contents{get(handles.popCSAlgo_2,'Value')};
fReplaceText([handles.currpath,filesep,'parameter_currRun_2.m'], sAlgo, 'cstype');
contents = cellstr(get(handles.popSparseTrafo_2,'String'));
sTrafo = contents{get(handles.popSparseTrafo_2,'Value')};
fReplaceText([handles.currpath,filesep,'parameter_currRun_2.m'], sTrafo, 'transformation');
edit([handles.currpath,filesep,'parameter_currRun_2.m']);


function [handles,imgScale] = fZeroRecon(handles)
% zero-padded reconstruction
if(~handles.measPara.prospectiveSub) % retrospective datasets
    if(~isempty(handles.dMask)) % subsample it first
        sDatadim = handles.measPara.dimension;
        iSparseDim = get(handles.edSparse,'Value');
        measPara = handles.measPara;
        if(strcmp(sDatadim,'2D'))
            switch iSparseDim
                case {1,4} % yPhase, zPhase | yPhase
                    iSparsePerm = [1, 2];
                    measPara.dimension = '2D';
                case 2 % yPhase, xFreq
                    iSparsePerm = [1, 3, 2];
                    measPara.dimension = '3D'; % faked workaround
                case {3,5} % xFreq, zPhase | xFreq
                    iSparsePerm = [2, 1];
                    measPara.dimension = '2D';
                case 6 % zPhase
                    errordlg('Chosen subsampling direction not possible for current dataset');
                    return;
            end 
        elseif(strcmp(sDatadim,'3D'))
            switch iSparseDim
                case 1 % yPhase, zPhase
                    iSparsePerm = [1, 2, 3];
                    measPara.dimension = '3D';
                case 2 % yPhase, xFreq
                    iSparsePerm = [1, 3, 2];
                    measPara.dimension = '3D';
                case 3 % xFreq, zPhase
                    iSparsePerm = [2, 1, 3];
                    measPara.dimension = '3D';
                case 4 % yPhase
                    iSparsePerm = [1, 2, 3];
                    measPara.dimension = '3D';
                case 5 % xFreq
                    iSparsePerm = [2, 1, 3];
                    measPara.dimension = '3D';
                case 6 % zPhase
                    iSparsePerm = [3, 1, 2];
                    measPara.dimension = '3D';
            end
        end
        kSpaceIn = cellfun(@(x) permute(x,iSparsePerm), handles.dKSpace, 'UniformOutput', false);
        if(strcmp(sDatadim,'3D') && strcmp(measPara.dimension,'2D')) % 3D -> 2D conversion was done
            kSpaceIn = cell2mat(shiftdim(kSpaceIn,-2));
            kSpaceIn = shiftdim(mat2cell(kSpaceIn, size(kSpaceIn,1), size(kSpaceIn,2), ones(1,size(kSpaceIn,3)), ones(1,size(kSpaceIn,4))),2);
        end
        iIndSparseDim = get(handles.edSparse,'Value');
        switch iIndSparseDim
            case 1 % yPhase, zPhase
                iPermMask = [1 3 2];
                iRepMask = [1 1 size(kSpaceIn{1},2)];
            case 2 % yPhase, xFreq
                iPermMask = [1 3 2];
                iRepMask = [1 1 size(kSpaceIn{1},2)];
            case 3 % xFreq, zPhase
                iPermMask = [1 3 2];
                iRepMask = [1 1 size(kSpaceIn{1},2)];
            case 4 % yPhase
                iPermMask = [1 2 3];
                iRepMask = [1 1 size(kSpaceIn{1},3)];
            case 5 % xFreq
                iPermMask = [2 1 3];
                iRepMask = [1 1 size(kSpaceIn{1},3)];
            case 6 % zPhase
                iPermMask = [1 3 2];
                iRepMask = [1 1 size(kSpaceIn{1},3)];        
        end

        lMask = repmat(handles.dMask,iRepMask);
        lMask = permute(lMask,iPermMask);
%         if(handles.dimMask(2) == 1) % 1D mask
%             lMask = repmat(handles.dMask, 1, 1, size(kSpaceIn{1},3));
%         else % 2D mask (just for 3D datasets, and faked 2D y-x)
%             lMask = repmat(handles.dMask, 1, 1, size(kSpaceIn{1},2));
%             lMask = permute(lMask, [1 3 2]);
%         end
        kSpaceIn = cellfun(@(x) x.*lMask, kSpaceIn, 'UniformOutput', false);

        para.cstype = 'Zero';
        para.transformation = 'fft';
        para.measPara = measPara;
        para.postproc.FreqOversamplingCorr = false;
        para.postproc.turnImage = false;
        para.prop.flagPlot = false;
        imgOut = CS_reconstruction(kSpaceIn, para);
        if(isequal(size(kSpaceIn{1}),size(imgOut)))
            imgOut = ipermute(imgOut,iSparsePerm);
        end
        if(strcmp(measPara.dimension,'2D'))
            imgScale = zeros(size(imgOut));
            for iI=1:size(imgScale,3)
                imgScale(:,:,iI) = scaleImg(imgOut(:,:,iI), handles.dRange);
            end
        else
            imgScale = scaleImg(imgOut, handles.dRange);
        end
    else
        imgScale = zeros(size(handles.dImg));
    end
else
    para.cstype = 'Zero';
    para.transformation = 'fft';
    para.measPara = handles.measPara;
    para.postproc.FreqOversamplingCorr = false;
    para.postproc.turnImage = false;
    para.prop.flagPlot = false;
    imgOut = CS_reconstruction(handles.dKSpace, para);
    if(strcmp(handles.measPara.dimension,'2D'))
        imgScale = zeros(size(imgOut));
        for iI=1:size(imgScale,3)
            imgScale(:,:,iI) = scaleImg(imgOut(:,:,iI), handles.dRange);
        end
    else
        imgScale = scaleImg(imgOut, handles.dRange);
    end
end


function popShowDiff_1_Callback(hObject, eventdata, handles)
% show difference images

iInd = get(hObject,'Value');
switch iInd
    case 1 % reference
        if(~handles.measPara.prospectiveSub) % retrospective datasets
            handles.dDiff_1 = scaleImg(handles.dImg, handles.dRange);
        else
            handles.dDiff_1 = zeros(size(handles.dImg));
        end
    case 2 % zero-padded
        [handles, dZero] = fZeroRecon(handles);
        handles.dDiff_1 = dZero;
    case 3 % zero-padded - recon_1
        [handles, dZero] = fZeroRecon(handles);
        handles.dDiff_1 = abs(dZero - handles.dOut_1);
    case 4 % reference - recon_1
        if(~handles.measPara.prospectiveSub) % retrospective datasets
            dRef = scaleImg(handles.dImg, handles.dRange);
        else
            dRef = zeros(size(handles.dImg));
        end
        handles.dDiff_1 = abs(dRef - handles.dOut_1);
    case 5 % recon_1 - recon_2
        handles.dDiff_1 = abs(handles.dOut_1 - handles.dOut_2);        
end
handles.hPlot.dDiff_1 = fPlotAxes(handles.dDiff_1, handles.iCurrSlice, handles.hPlot.dDiff_1, [], handles.dRange);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popShowDiff_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popShowDiff_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popShowDiff_2.
function popShowDiff_2_Callback(hObject, eventdata, handles)
% show difference images

iInd = get(hObject,'Value');
switch iInd
    case 1 % reference
        if(~handles.measPara.prospectiveSub) % retrospective datasets
            handles.dDiff_2 = scaleImg(handles.dImg, handles.dRange);
        else
            handles.dDiff_2 = zeros(size(handles.dImg));
        end
    case 2 % zero-padded
        [handles, dZero] = fZeroRecon(handles);
        handles.dDiff_2 = dZero;
    case 3 % zero-padded - recon_1
        [handles, dZero] = fZeroRecon(handles);
        handles.dDiff_2 = abs(dZero - handles.dOut_2);
    case 4 % reference - recon_2
        if(~handles.measPara.prospectiveSub) % retrospective datasets
            dRef = scaleImg(handles.dImg, handles.dRange);
        else
            dRef = zeros(size(handles.dImg));
        end
        handles.dDiff_2 = abs(dRef - handles.dOut_2);
    case 5 % recon_1 - recon_2
        handles.dDiff_2 = abs(handles.dOut_1 - handles.dOut_2);        
end
handles.hPlot.dDiff_2 = fPlotAxes(handles.dDiff_2, handles.iCurrSlice, handles.hPlot.dDiff_2, [], handles.dRange);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popShowDiff_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popShowDiff_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function hImshow = fPlotAxes(dData, iCurrSlice, hImshow, hAxes, dRange)
% general image plot function

% handles.(handles.idnames{2,iID}) = scaleImg(dData, handles.dRange); % 2D!
% axes(handles.(handles.idnames{1,iID}));
% set(handles.hPlot.(handles.idnames{2,iID}), 'CData', handles.(handles.idnames{2,iID}));
% colormap(handles.(handles.idnames{1,iID}),'gray');
if(iCurrSlice > size(dData,3))
    iCurrSlice = size(dData,3);
elseif(iCurrSlice < 1)
    iCurrSlice = 1;
end
if(isempty(hAxes)) % just update
    set(hImshow, 'CData', dData(:,:,iCurrSlice));
else % refresh axes
    cla(hAxes);
    axes(hAxes);
    hImshow = imshow(dData(:,:,iCurrSlice), dRange);
    colormap(hAxes,'gray');
end


function [kSpace, kSpaceShow, dImg, mask, measPara, sInfoTxt] = fLoadDataset(sPath, sparseDim)
% general load data function
if(~exist(sPath,'file'))
    return;
end
kSpace = [];
mask = [];
measPara = [];
load(sPath);
if(isempty(kSpace))
    return;
end
if(strcmp(measPara.dimension, '2D'))
    fftDims = 1:2;
else
    fftDims = 1:3;
end
if(isempty(measPara))
    measPara.dim = [size(kSpace{1,1,1,1},1), size(kSpace{1,1,1,1},2), size(kSpace{1,1,1,1},3), size(kSpace{1,1,1,1},4), size(kSpace,2)];
    switch ndims(kSpace{1,1,1,1})
        case 2 % 2D
            measPara.dimension = '2D';
            measPara.dim(3) = size(kSpace,1);
            fftDims = 1:2;
        case 3 % 3D
            measPara.dimension = '3D';
            fftDims = 1:3;
        case 4 % 4D
            if(measPara.dim(3) == 1) % 2Dt
                measPara.dimension = '2D';
                kSpace = cellfun(@squeeze, kSpace, 'UniformOutput', false);
                fftDims = 1:2;
            else % 4D
                measPara.dimension = '4D';
                fftDims = 1:3;
            end
    end
end
sInfoTxt{1} = sprintf('size (yPhase / xFreq / zPhase|slice / time / channel): %d / %d / %d / %d / %d', measPara.dim(1), measPara.dim(2), measPara.dim(3), measPara.dim(4), measPara.dim(5));
switch sparseDim
    case 1 % yPhase, zPhase
        sInfoTxt{2} = sprintf('mask dim: %d x %d', measPara.dim(1), measPara.dim(3));
    case 2 % yPhase, xFreq
        sInfoTxt{2} = sprintf('mask dim: %d x %d', measPara.dim(1), measPara.dim(2));
    case 3 % xFreq, zPhase
        sInfoTxt{2} = sprintf('mask dim: %d x %d', measPara.dim(2), measPara.dim(3));
    case 4 % yPhase
        sInfoTxt{2} = sprintf('mask dim: %d x 1', measPara.dim(1));
    case 5 % xFreq
        sInfoTxt{2} = sprintf('mask dim: %d x 1', measPara.dim(2));
    case 6 % zPhase
        sInfoTxt{2} = sprintf('mask dim: %d x 1', measPara.dim(3));
end
kSpaceShow = cell2mat(shiftdim(kSpace,-2));
dImg = squeeze(sqrt(sum(abs(ifftnshift(kSpaceShow,fftDims)).^2,4)));
kSpaceShow = squeeze(abs(sum(kSpaceShow,4)));

function scaledImg = scaleImg(img,range)
scaledImg = ((img - min(img(:))) * (range(2)-range(1)))./(max(img(:)) - min(img(:)));

function fReplaceText(fFilename, sText, sWhatToLookFor)
%     if(nargin < 3)
%         sWhatReplaced = 'InitialTransformParametersFileName';
%     end
fid = fopen(fFilename,'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
    if(ischar(A{i}) && ~isempty(regexp(A{i},[sprintf('^%s',sWhatToLookFor),'\w*'],'once')))
        if(ischar(sText))
            A{i} = sprintf('%s = ''%s'';',sWhatToLookFor,sText);
        else % double array 
            sTmp = '[';
            for iJ=1:size(sText,1)
                for iK=1:size(sText,2)
                    sTmp = [sTmp, sprintf('%d ',sText(iJ,iK))];
                end
                if(iJ == size(sText,1))
                    sTmp = [sTmp, ']'];
                else
                    sTmp = [sTmp, ';'];
                end
            end
            A{i} = sprintf('%s = %s;',sWhatToLookFor,sTmp);
        end
    end
end
fclose(fid);
fid = fopen(fFilename, 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break;
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);

function [ samplingMask ] = readSamplingMaskFromFile(sFilename)
%This function reads out a sampling pattern from a .txt-file

if(nargin < 1)
    file = fopen('samplingPattern.txt','r');
else
    file = fopen(sFilename,'r');
end

endOfDoc = false;
i = 0;

while(~endOfDoc)
    actualLine = fgetl(file);
    
    if(actualLine == -1)
        endOfDoc = true;
        break;
    end
    i = i + 1;
    lines{i} = actualLine;
    
    for j=1:numel(lines{i})
        actualChar = lines{i}(j);

        if(actualChar == '-')
            samplingMask(i, j) = 0;
        elseif (actualChar == 'x' || actualChar == '*' || actualChar == 'O' || actualChar == '0')
            samplingMask(i, j) = 1;
        else
            error('not working, because there is an unknown charakter');
        end
    end
end

if(exist('drawMask','file'))    
    drawMask(samplingMask);
end

fclose(file);


% --- Executes on scroll wheel click while the figure is in focus.
function CS_LAB_GUI_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to CS_LAB_GUI (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
% get scroll events

lMoving = true;
if eventdata.VerticalScrollCount < 0
    handles.iCurrSlice = max([1 handles.iCurrSlice - 1]);
elseif(eventdata.VerticalScrollCount > 0)
    handles.iCurrSlice = min([size(handles.dInput, 3) handles.iCurrSlice + 1]);
else 
    lMoving = false;
end

if(lMoving)
    for iID=1:size(handles.idnames,2)
        handles.hPlot.(handles.idnames{2,iID}) = fPlotAxes(handles.(handles.idnames{2,iID}), handles.iCurrSlice, handles.hPlot.(handles.idnames{2,iID}), [], []);
    end
end

guidata(hObject, handles);


function pbRun_1_Callback(hObject, eventdata, handles)
% execute algorithm1

if(nnz(handles.dRecon_1) > 0)
    choice = questdlg('Existing reconstruction will be deleted. Do you want to continue?', 'Deleting results', 'Yes', 'No', 'Yes');
    switch choice
        case 'Yes'
            % NOP
        case 'No'
            return;
    end
end

measPara = handles.measPara;
if(handles.iEsp(1) == 1)
    espresso.state = false;
    espresso.direction = 'off';
else
    espresso.state = true;
    espresso.direction = 'on'; % workaround, actual direction info not needed for recon (will be extracted automatically)
end
espresso.pfn = handles.iEsp(1);
% get sparse dimensions
sDatadim = handles.measPara.dimension;
iSparseDim = get(handles.edSparse,'Value');
if(iSparseDim <= 3)
    iSparseDimensionality = 2; % 2D sparse trafo necessary
else
    iSparseDimensionality = 1; % just 1D sparse trafo necessary
end
iDatadim = str2num(sDatadim(1));
set(handles.popTrafoDim_1,'Value', max([min([get(handles.popTrafoDim_1,'Value'),iDatadim]), iSparseDimensionality])); % ensure correct dimensionality
iTrafoDim = get(handles.popTrafoDim_1,'Value');
sTrafodim = '[1 2 1 0; '; % ATTENTION: faked workaround (so no adaptation inside CS_recon is necessary) trafoDim does no longer reflect [y x z t] due to permutation of iSparsePerm => switch sparse dimensions into X-marked positions trafoDim = [X 0 X 0]
if(strcmp(sDatadim,'2D'))
    switch iSparseDim
        case {1,4} % yPhase, zPhase | yPhase => zPhase = 1, same as 4
            iSparsePerm = [1, 2];
            measPara.dimension = '2D';
            if(iTrafoDim == 1)
                sTrafodim = [sTrafodim, '1 0 0 0];'];
            elseif(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 1 0 0];'];
            end              
        case 2 % yPhase, xFreq
            iSparsePerm = [1, 3, 2];
            measPara.dimension = '3D'; % faked workaround
            sTrafodim = '[1 0 1 0; 1 0 1 0];';
%             sTrafodim = [sTrafodim, '1 0 1 0];'];
        case {3,5} % xFreq, zPhase | xFreq => zPhase = 1, same as 5
            iSparsePerm = [2, 1]; 
            measPara.dimension = '2D';
            if(iTrafoDim == 1)
                sTrafodim = [sTrafodim, '1 0 0 0];'];
            elseif(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 1 0 0];'];
            end
        case 6 % zPhase
            errordlg('Chosen subsampling direction not possible for current dataset');
            return;
    end 
elseif(strcmp(sDatadim,'3D'))
    switch iSparseDim
        case 1 % yPhase, zPhase
            iSparsePerm = [1, 2, 3];
            measPara.dimension = '3D';
            if(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 0 1 0];'];
            elseif(iTrafoDim == 3)
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
        case 2 % yPhase, xFreq
            iSparsePerm = [1, 3, 2];
            measPara.dimension = '3D';
            if(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 0 1 0];'];
            elseif(iTrafoDim == 3)
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
        case 3 % xFreq, zPhase
            iSparsePerm = [2, 1, 3];
            measPara.dimension = '3D';
            if(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 0 1 0];'];
            elseif(iTrafoDim == 3)
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
        case 4 % yPhase
            iSparsePerm = [1, 2, 3];
            measPara.dimension = '2D';
            if(iTrafoDim == 1)
                sTrafodim = [sTrafodim, '1 0 0 0];'];
            elseif(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 1 0 0];'];
            elseif(iTrafoDim == 3)
                measPara.dimension = '3D';
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
        case 5 % xFreq
            iSparsePerm = [2, 1, 3];
            measPara.dimension = '2D';
            if(iTrafoDim == 1)
                sTrafodim = [sTrafodim, '1 0 0 0];'];
            elseif(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 1 0 0];'];
            elseif(iTrafoDim == 3)
                measPara.dimension = '3D';
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
        case 6 % zPhase
            iSparsePerm = [3, 1, 2];
            measPara.dimension = '2D';
            if(iTrafoDim == 1)
                sTrafodim = [sTrafodim, '1 0 0 0];'];
            elseif(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 1 0 0];'];
            elseif(iTrafoDim == 3)
                measPara.dimension = '3D';
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
    end
end
kSpaceIn = cellfun(@(x) permute(x,iSparsePerm), handles.dKSpace, 'UniformOutput', false);
if(strcmp(sDatadim,'3D') && strcmp(measPara.dimension,'2D')) % 3D -> 2D conversion was done    
    kSpaceIn = cell2mat(shiftdim(kSpaceIn,-2));
    kSpaceIn = ifftnshift(kSpaceIn,3); % zPhase -> slice encoding
    kSpaceIn = shiftdim(mat2cell(kSpaceIn, size(kSpaceIn,1), size(kSpaceIn,2), ones(1,size(kSpaceIn,3)), ones(1,size(kSpaceIn,4))),2);
end
% prepare k-space input
if(~measPara.prospectiveSub) % retrospective datasets
    if(isempty(handles.dMask)) % generate subsampling mask first
        handles = fCreateMask(handles);
    end
    iIndSparseDim = get(handles.edSparse,'Value');
    switch iIndSparseDim
        case 1 % yPhase, zPhase
            iPermMask = [1 3 2];
            iRepMask = [1 1 size(kSpaceIn{1},2)];
            iSparseSize = [handles.measPara.dim(1), handles.measPara.dim(3)];
        case 2 % yPhase, xFreq
            iPermMask = [1 3 2];
            iRepMask = [1 1 size(kSpaceIn{1},2)];
            iSparseSize = [handles.measPara.dim(1), handles.measPara.dim(2)];
        case 3 % xFreq, zPhase
            iPermMask = [1 3 2];
            iRepMask = [1 1 size(kSpaceIn{1},2)];
            iSparseSize = [handles.measPara.dim(2), handles.measPara.dim(3)];
        case 4 % yPhase
            iPermMask = [1 2 3];
            iRepMask = [1 1 size(kSpaceIn{1},3)];
            iSparseSize = [handles.measPara.dim(1), 1];
        case 5 % xFreq
            iPermMask = [2 1 3];
            iRepMask = [1 1 size(kSpaceIn{1},3)];
            iSparseSize = [handles.measPara.dim(2), 1];
        case 6 % zPhase
            iPermMask = [1 2 3];
            iRepMask = [1 1 size(kSpaceIn{1},3)];
            iSparseSize = [handles.measPara.dim(3), 1];
    end
    
    lMask = handles.dMask;
    if(ndims(lMask) ~= length(iSparseSize) || any(size(lMask) ~= iSparseSize))
        lMask = lMask(1:iSparseSize(1),1:iSparseSize(2));
    end
    lMask = repmat(lMask,iRepMask);
    lMask = permute(lMask,iPermMask);
    
%     if(handles.dimMask(2) == 1) % 1D mask
%         lMask = repmat(handles.dMask, 1, 1, size(kSpaceIn{1},3)); % handles.dMask already repmatted towards other shown image dim
%     end
%     else % 2D mask (just for 3D datasets, and faked 2D y-x)
%         lMask = repmat(handles.dMask, 1, 1, size(kSpaceIn{1},2));
%         lMask = permute(lMask, [1 3 2]);
%     end        
    kSpaceIn = cellfun(@(x) x.*lMask, kSpaceIn, 'UniformOutput', false);  
    if(nnz(handles.dDiff_1) == 0)
        handles.dDiff_1 = scaleImg(handles.dImg,handles.dRange);
    end
end

% create parameter file
if(~exist([handles.currpath,filesep,'parameter_currRun_1.m'],'file'))
    copyfile([handles.currpath,filesep,'parameters_default_GUI.m'], [handles.currpath,filesep,'parameter_currRun_1.m']);
end
contents = cellstr(get(handles.popCSAlgo_1,'String'));
sAlgo = contents{get(handles.popCSAlgo_1,'Value')};
fReplaceText([handles.currpath,filesep,'parameter_currRun_1.m'], sAlgo, 'cstype');
contents = cellstr(get(handles.popSparseTrafo_1,'String'));
sTrafo = contents{get(handles.popSparseTrafo_1,'Value')};
fReplaceText([handles.currpath,filesep,'parameter_currRun_1.m'], sTrafo, 'transformation');
fReplaceText([handles.currpath,filesep,'parameter_currRun_1.m'], str2num(sTrafodim), 'trafo.trafodim');
% add measPara struct
fid = fopen([handles.currpath,filesep,'parameter_currRun_1.m'],'a+');
fprintf(fid,'\n');
cFieldnames = fieldnames(measPara);
for iI=1:length(cFieldnames)
    switch cFieldnames{iI}
        case 'dimension'
            sFormat = 'measPara.%s = ''%s'';\n';
        case 'dim'
            sFormat = 'measPara.%s = [%d %d %d %d %d];\n';
        case 'LCall'
            sFormat = 'measPara.%s = [%d %d %d %d];\n';
        otherwise 
            continue;
    end
    fprintf(fid, sFormat, cFieldnames{iI}, measPara.(cFieldnames{iI}));
end
% add espresso info
cFieldnames = fieldnames(espresso);
sTmp = {'false','true'};
for iI=1:length(cFieldnames)
    switch cFieldnames{iI}
        case 'pfn'
            sFormat = 'espresso.%s = %f;\n';
        case 'state'            
            sFormat = 'espresso.%s = logical(%d);\n';
        case 'direction'
            sFormat = 'espresso.%s = ''%s'';\n';
    end
    fprintf(fid, sFormat, cFieldnames{iI}, espresso.(cFieldnames{iI}));
end
fclose(fid);

handles.dRecon_1 = CS_reconstruction(kSpaceIn, [handles.currpath,filesep,'parameter_currRun_1.m']); 
cd(handles.currpath);
if(strcmp(measPara.dimension,'2D'))
    handles.dOut_1 = zeros(size(handles.dRecon_1));
    for iI=1:size(handles.dRecon_1,3)
        handles.dOut_1(:,:,iI) = scaleImg(handles.dRecon_1(:,:,iI), handles.dRange);
    end
else
    handles.dOut_1 = scaleImg(handles.dRecon_1, handles.dRange);
end
if(isequal(size(kSpaceIn{1}),size(handles.dRecon_1)))
    handles.dRecon_1 = ipermute(handles.dRecon_1,iSparsePerm);
    handles.dOut_1 = ipermute(handles.dOut_1,iSparsePerm);
end
handles.hPlot.dOut_1 = fPlotAxes(handles.dOut_1, handles.iCurrSlice, handles.hPlot.dOut_1, [], handles.dRange);
handles.hPlot.dDiff_1 = fPlotAxes(handles.dDiff_1, handles.iCurrSlice, handles.hPlot.dDiff_1, [], handles.dRange);
cd(handles.currpath);

guidata(hObject, handles);


function pbRun_2_Callback(hObject, eventdata, handles)
% execute algorithm2

if(nnz(handles.dRecon_2) > 0)
    choice = questdlg('Existing reconstruction will be deleted. Do you want to continue?', 'Deleting results', 'Yes', 'No', 'Yes');
    switch choice
        case 'Yes'
            % NOP
        case 'No'
            return;
    end
end

measPara = handles.measPara;
if(handles.iEsp(1) == 1)
    espresso.state = false;
    espresso.direction = 'off';
else
    espresso.state = true;
    espresso.direction = 'on'; % workaround, actual direction info not needed for recon (will be extracted automatically)
end
espresso.pfn = handles.iEsp(1);
% get sparse dimensions
sDatadim = handles.measPara.dimension;
iSparseDim = get(handles.edSparse,'Value');
if(iSparseDim <= 3)
    iSparseDimensionality = 2; % 2D sparse trafo necessary
else
    iSparseDimensionality = 1; % just 1D sparse trafo necessary
end
iDatadim = str2num(sDatadim(1));
set(handles.popTrafoDim_2,'Value', max([min([get(handles.popTrafoDim_2,'Value'),iDatadim]), iSparseDimensionality]));
iTrafoDim = get(handles.popTrafoDim_2,'Value');
sTrafodim = '[1 2 1 0; '; % ATTENTION: faked workaround (so no adaptation inside CS_recon is necessary) trafoDim does no longer reflect [y x z t] due to permutation of iSparsePerm => switch sparse dimensions into X-marked positions trafoDim = [X 0 X 0]
if(strcmp(sDatadim,'2D'))
    switch iSparseDim
        case {1,4} % yPhase, zPhase | yPhase => zPhase = 1, same as 4
            iSparsePerm = [1, 2];
            measPara.dimension = '2D';
            if(iTrafoDim == 1)
                sTrafodim = [sTrafodim, '1 0 0 0];'];
            elseif(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 1 0 0];'];
            end              
        case 2 % yPhase, xFreq
            iSparsePerm = [1, 3, 2];
            measPara.dimension = '3D'; % faked workaround
            sTrafodim = [sTrafodim, '1 0 1 0];'];
        case {3,5} % xFreq, zPhase | xFreq => zPhase = 1, same as 5
            iSparsePerm = [2, 1]; 
            measPara.dimension = '2D';
            if(iTrafoDim == 1)
                sTrafodim = [sTrafodim, '1 0 0 0];'];
            elseif(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 1 0 0];'];
            end
        case 6 % zPhase
            errordlg('Chosen subsampling direction not possible for current dataset');
            return;
    end 
elseif(strcmp(sDatadim,'3D'))
    switch iSparseDim
        case 1 % yPhase, zPhase
            iSparsePerm = [1, 2, 3];
            measPara.dimension = '3D';
            if(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 0 1 0];'];
            elseif(iTrafoDim == 3)
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
        case 2 % yPhase, xFreq
            iSparsePerm = [1, 3, 2];
            measPara.dimension = '3D';
            if(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 0 1 0];'];
            elseif(iTrafoDim == 3)
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
        case 3 % xFreq, zPhase
            iSparsePerm = [2, 1, 3];
            measPara.dimension = '3D';
            if(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 0 1 0];'];
            elseif(iTrafoDim == 3)
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
        case 4 % yPhase
            iSparsePerm = [1, 2, 3];
            measPara.dimension = '2D';
            if(iTrafoDim == 1)
                sTrafodim = [sTrafodim, '1 0 0 0];'];
            elseif(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 1 0 0];'];
            elseif(iTrafoDim == 3)
                measPara.dimension = '3D';
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
        case 5 % xFreq
            iSparsePerm = [2, 1, 3];
            measPara.dimension = '2D';
            if(iTrafoDim == 1)
                sTrafodim = [sTrafodim, '1 0 0 0];'];
            elseif(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 1 0 0];'];
            elseif(iTrafoDim == 3)
                measPara.dimension = '3D';
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
        case 6 % zPhase
            iSparsePerm = [3, 1, 2];
            measPara.dimension = '2D';
            if(iTrafoDim == 1)
                sTrafodim = [sTrafodim, '1 0 0 0];'];
            elseif(iTrafoDim == 2)
                sTrafodim = [sTrafodim, '1 1 0 0];'];
            elseif(iTrafoDim == 3)
                measPara.dimension = '3D';
                sTrafodim = [sTrafodim, '1 1 1 0];'];
            end
    end
end
kSpaceIn = cellfun(@(x) permute(x,iSparsePerm), handles.dKSpace, 'UniformOutput', false);
if(strcmp(sDatadim,'3D') && strcmp(measPara.dimension,'2D')) % 3D -> 2D conversion was done
    kSpaceIn = cell2mat(shiftdim(kSpaceIn,-2));
    kSpaceIn = ifftnshift(kSpaceIn,3); % zPhase -> slice encoding
    kSpaceIn = shiftdim(mat2cell(kSpaceIn, size(kSpaceIn,1), size(kSpaceIn,2), ones(1,size(kSpaceIn,3)), ones(1,size(kSpaceIn,4))),2);
end
% prepare k-space input
if(~measPara.prospectiveSub) % retrospective datasets
    if(isempty(handles.dMask)) % generate subsampling mask first
        handles = fCreateMask(handles);
    end
    iIndSparseDim = get(handles.edSparse,'Value');
    switch iIndSparseDim
        case 1 % yPhase, zPhase
            iPermMask = [1 3 2];
            iRepMask = [1 1 size(kSpaceIn{1},2)];
            iSparseSize = [handles.measPara.dim(1), handles.measPara.dim(3)];
        case 2 % yPhase, xFreq
            iPermMask = [1 3 2];
            iRepMask = [1 1 size(kSpaceIn{1},2)];
            iSparseSize = [handles.measPara.dim(1), handles.measPara.dim(2)];
        case 3 % xFreq, zPhase
            iPermMask = [1 3 2];
            iRepMask = [1 1 size(kSpaceIn{1},2)];
            iSparseSize = [handles.measPara.dim(2), handles.measPara.dim(3)];
        case 4 % yPhase
            iPermMask = [1 2 3];
            iRepMask = [1 1 size(kSpaceIn{1},3)];
            iSparseSize = [handles.measPara.dim(1), 1];
        case 5 % xFreq
            iPermMask = [2 1 3];
            iRepMask = [1 1 size(kSpaceIn{1},3)];
            iSparseSize = [handles.measPara.dim(2), 1];
        case 6 % zPhase
            iPermMask = [1 2 3];
            iRepMask = [1 1 size(kSpaceIn{1},3)];    
            iSparseSize = [handles.measPara.dim(3), 1];
    end
    
    lMask = handles.dMask;
    if(ndims(lMask) ~= length(iSparseSize) || any(size(lMask) ~= iSparseSize))
        lMask = lMask(1:iSparseSize(1),1:iSparseSize(2));
    end
    lMask = repmat(lMask,iRepMask);
    lMask = permute(lMask,iPermMask);
    
%     if(handles.dimMask(2) == 1) % 1D mask
%         lMask = repmat(handles.dMask, 1, size(kSpaceIn{1},2), size(kSpaceIn{1},3));
%     else % 2D mask (just for 3D datasets, and faked 2D y-x)
%         lMask = repmat(handles.dMask, 1, 1, size(kSpaceIn{1},2));
%         lMask = permute(lMask, [1 3 2]);
%     end        
    kSpaceIn = cellfun(@(x) x.*lMask, kSpaceIn, 'UniformOutput', false);    
    if(nnz(handles.dDiff_2) == 0)
        handles.dDiff_2 = scaleImg(handles.dImg,handles.dRange);
    end
end

% create parameter file
if(~exist([handles.currpath,filesep,'parameter_currRun_2.m'],'file'))
    copyfile([handles.currpath,filesep,'parameters_default_GUI.m'], [handles.currpath,filesep,'parameter_currRun_2.m']);
end
contents = cellstr(get(handles.popCSAlgo_2,'String'));
sAlgo = contents{get(handles.popCSAlgo_2,'Value')};
fReplaceText([handles.currpath,filesep,'parameter_currRun_2.m'], sAlgo, 'cstype');
contents = cellstr(get(handles.popSparseTrafo_2,'String'));
sTrafo = contents{get(handles.popSparseTrafo_2,'Value')};
fReplaceText([handles.currpath,filesep,'parameter_currRun_2.m'], sTrafo, 'transformation');
fReplaceText([handles.currpath,filesep,'parameter_currRun_2.m'], str2num(sTrafodim), 'trafo.trafodim');
% add measPara struct
fid = fopen([handles.currpath,filesep,'parameter_currRun_2.m'],'a+');
fprintf(fid,'\n');
cFieldnames = fieldnames(measPara);
for iI=1:length(cFieldnames)
    switch cFieldnames{iI}
        case 'dimension'
            sFormat = 'measPara.%s = ''%s'';\n';
        case 'dim'
            sFormat = 'measPara.%s = [%d %d %d %d %d];\n';
        case 'LCall'
            sFormat = 'measPara.%s = [%d %d %d %d];\n';
        otherwise 
            continue;
    end
    fprintf(fid, sFormat, cFieldnames{iI}, measPara.(cFieldnames{iI}));
end
% add espresso info
cFieldnames = fieldnames(espresso);
sTmp = {'false','true'};
for iI=1:length(cFieldnames)
    switch cFieldnames{iI}
        case 'pfn'
            sFormat = 'espresso.%s = %f;\n';
        case 'state'            
            sFormat = 'espresso.%s = sTmp{%d+1};\n';
        case 'direction'
            sFormat = 'espresso.%s = %s;\n';
    end
    fprintf(fid, sFormat, cFieldnames{iI}, espresso.(cFieldnames{iI}));
end
fclose(fid);

handles.dRecon_2 = CS_reconstruction(kSpaceIn, [handles.currpath,filesep,'parameter_currRun_2.m']);
cd(handles.currpath);
if(strcmp(measPara.dimension,'2D'))
    handles.dOut_2 = zeros(size(handles.dRecon_2));
    for iI=1:size(handles.dRecon_2,3)
        handles.dOut_2(:,:,iI) = scaleImg(handles.dRecon_2(:,:,iI), handles.dRange);
    end
else
    handles.dOut_2 = scaleImg(handles.dRecon_2, handles.dRange);
end
if(isequal(size(kSpaceIn{1}),size(handles.dRecon_2)))
    handles.dRecon_2 = ipermute(handles.dRecon_2,iSparsePerm);
    handles.dOut_2 = ipermute(handles.dOut_2,iSparsePerm);
end
handles.hPlot.dOut_2 = fPlotAxes(handles.dOut_2, handles.iCurrSlice, handles.hPlot.dOut_2, [], handles.dRange);
handles.hPlot.dDiff_2 = fPlotAxes(handles.dDiff_2, handles.iCurrSlice, handles.hPlot.dDiff_2, [], handles.dRange);
cd(handles.currpath);
guidata(hObject, handles);


function edSparse_Callback(hObject, eventdata, handles)
% change sparse dimensions for subsampling mask
switch get(hObject,'Value')
    case 1 % yPhase, zPhase
        sInfoTxt = sprintf('mask dim: %d x %d', handles.measPara.dim(1), handles.measPara.dim(3));
    case 2 % yPhase, xFreq
        sInfoTxt = sprintf('mask dim: %d x %d', handles.measPara.dim(1), handles.measPara.dim(2));
    case 3 % xFreq, zPhase
        sInfoTxt = sprintf('mask dim: %d x %d', handles.measPara.dim(2), handles.measPara.dim(3));
    case 4 % yPhase
        sInfoTxt = sprintf('mask dim: %d x 1', handles.measPara.dim(1));
    case 5 % xFreq
        sInfoTxt = sprintf('mask dim: %d x 1', handles.measPara.dim(2));
    case 6 % zPhase
        sInfoTxt = sprintf('mask dim: %d x 1', handles.measPara.dim(3));
end
set(handles.txtMask,'String',sInfoTxt);


% --- Executes during object creation, after setting all properties.
function edSparse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edSparse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pbSave_1_Callback(hObject, eventdata, handles)
% save results of algo1
[FileName,PathName] = uiputfile('*.mat','Save result file');
if(isempty(FileName))
    return;
end
dRecon = handles.dRecon_1;
dOut = handles.dOut_1;
dMask = handles.dMask;
% save reference as dDiff_1
iIndDataset = get(handles.popDatasets,'Value');
sInputFile = handles.cDatasets{2,iIndDataset};
if(~handles.measPara.prospectiveSub) % retrospective
    dDiff = handles.dImg;
else
    dDiff = zeros(size(handles.dImg));
end
save([PathName,filesep,FileName], 'sInputFile', 'dRecon', 'dOut', 'dMask', 'dDiff');


function pbLoad_1_Callback(hObject, eventdata, handles)
% load results of algo1
[FileName,PathName] = uigetfile('*.mat','Select result file');
if(isempty(FileName))
    return;
end
load([PathName,filesep,FileName]);
if(~exist('dRecon','var'))
    return;
end
handles.dRecon_1 = dRecon;
% handles.dOut_1 = scaleImg(handles.dRecon_1, handles.dRange);
handles.dOut_1 = dOut;
handles.dMask = dMask;
handles.dDiff_1 = dDiff;

[~,sFilename] = fileparts(sInputFile);
handles.cDatasets(:,end+1) = {sFilename; sInputFile};
iIndDataset = size(handles.cDatasets,2);
set(handles.popDatasets,'String', handles.cDatasets(1,:));
set(handles.popDatasets,'Value',iIndDataset);
[handles.dKSpace, handles.dKSpaceShow, handles.dImg, dMaskLoad, handles.measPara, sInfoTxt] = fLoadDataset(handles.cDatasets{2,iIndDataset}, get(handles.edSparse,'Value'));
if(~isempty(dMaskLoad)) % mask was supplied
    handles.measPara.prospectiveSub = true;
    handles.dimMask = size(handles.dMask);
    set(handles.edSparse,'Value',handles.measPara.iIndSparseDim);   
    set(handles.edSparse,'Enable','off');
else
    handles.measPara.prospectiveSub = false;
    set(handles.edSparse,'Enable','on');
end
if(size(handles.dImg,3) == 1) % just 2D image
    set(handles.popTrafoDim_1,'Value',2);
    set(handles.popTrafoDim_2,'Value',2);
else % 3D image
    set(handles.popTrafoDim_1,'Value',3);
    set(handles.popTrafoDim_2,'Value',3);
end
% get shown parameter
iInd = get(handles.popShowInput, 'Value');
switch iInd
    case 1 % image
        handles.dInput = scaleImg(handles.dImg, handles.dRange);        
    case 2 % kspace
        handles.dInput = scaleImg(handles.dKSpaceShow, handles.dRange);    
end
handles.iCurrSlice = round(size(handles.dInput,3)/2);
for iI=1:size(handles.idnames,2)
    handles.hPlot.(handles.idnames{2,iI}) = fPlotAxes(handles.(handles.idnames{2,iI}), handles.iCurrSlice, handles.hPlot.(handles.idnames{2,iI}), handles.(handles.idnames{1,iI}), handles.dRange);
end
set(handles.txtInfo,'String',sInfoTxt{1});
set(handles.txtMask,'String',sInfoTxt{2});
guidata(hObject, handles);        


function pbSave_2_Callback(hObject, eventdata, handles)
% save results of algo2
[FileName,PathName] = uiputfile('*.mat','Save result file');
if(isempty(FileName))
    return;
end
dRecon = handles.dRecon_2;
dOut = handles.dOut_2;
dMask = handles.dMask;
% save reference as dDiff_2
iIndDataset = get(handles.popDatasets,'Value');
sInputFile = handles.cDatasets{2,iIndDataset};
if(~handles.measPara.prospectiveSub) % retrospective
    dDiff = handles.dImg;
else
    dDiff = zeros(size(handles.dImg));
end
save([PathName,filesep,FileName], 'sInputFile', 'dRecon', 'dOut', 'dMask', 'dDiff');


function pbLoad_2_Callback(hObject, eventdata, handles)
% load results of algo2
[FileName,PathName] = uigetfile('*.mat','Select result file');
if(isempty(FileName))
    return;
end
load([PathName,filesep,FileName]);
if(~exist('dRecon','var'))
    return;
end
handles.dRecon_2 = dRecon;
% handles.dOut_2 = scaleImg(handles.dRecon_2, handles.dRange);
handles.dOut_2 = dOut;
handles.dMask = dMask;
handles.dDiff_2 = dDiff;

[~,sFilename] = fileparts(sInputFile);
handles.cDatasets(:,end+1) = {sFilename; sInputFile};
iIndDataset = size(handles.cDatasets,2);
set(handles.popDatasets,'String', handles.cDatasets(1,:));
set(handles.popDatasets,'Value',iIndDataset);
[handles.dKSpace, handles.dKSpaceShow, handles.dImg, dMaskLoad, handles.measPara, sInfoTxt] = fLoadDataset(handles.cDatasets{2,iIndDataset}, get(handles.edSparse,'Value'));
if(~isempty(dMaskLoad)) % mask was supplied
    handles.measPara.prospectiveSub = true;
    handles.dimMask = size(handles.dMask);
    set(handles.edSparse,'Value',handles.measPara.iIndSparseDim);   
    set(handles.edSparse,'Enable','off');
else
    handles.measPara.prospectiveSub = false;
    set(handles.edSparse,'Enable','on');
end
if(size(handles.dImg,3) == 1) % just 2D image
    set(handles.popTrafoDim_1,'Value',2);
    set(handles.popTrafoDim_2,'Value',2);
else % 3D image
    set(handles.popTrafoDim_1,'Value',3);
    set(handles.popTrafoDim_2,'Value',3);
end
% get shown parameter
iInd = get(handles.popShowInput, 'Value');
switch iInd
    case 1 % image
        handles.dInput = scaleImg(handles.dImg, handles.dRange);        
    case 2 % kspace
        handles.dInput = scaleImg(handles.dKSpaceShow, handles.dRange);    
end
handles.iCurrSlice = round(size(handles.dInput,3)/2);
for iI=1:size(handles.idnames,2)
    handles.hPlot.(handles.idnames{2,iI}) = fPlotAxes(handles.(handles.idnames{2,iI}), handles.iCurrSlice, handles.hPlot.(handles.idnames{2,iI}), handles.(handles.idnames{1,iI}), handles.dRange);
end
set(handles.txtInfo,'String',sInfoTxt{1});
set(handles.txtMask,'String',sInfoTxt{2});
guidata(hObject, handles);        


function CS_LAB_GUI_WindowButtonDownFcn(hObject, eventdata, handles)
% Save starting parameters for brightness/contrast scaling
handles.dPosStart = get(gca, 'CurrentPoint');

handles.FButtonDown = 1;
handles.colMin = handles.dRange(1);
handles.colMax = handles.dRange(2);
guidata(hObject, handles)


function CS_LAB_GUI_WindowButtonMotionFcn(hObject, eventdata, handles)
% apply brightness/contrast scaling
try
    if handles.FButtonDown
        switch get(hObject, 'SelectionType')
            case 'extend'
                iD = get(gca, 'CurrentPoint') - handles.dPosStart; % Mouse distance travelled since button down
                
                % contrast and brightness
                handles.colWidth  = handles.colMax-handles.colMin;
                handles.colWidth  = handles.colMax.*exp(-iD(1,2)*0.02);
                handles.colCenter = (handles.colMax+handles.colMin)/2;
                handles.colCenter = handles.colCenter.*exp(iD(1,1)*0.02);
                handles.dRange  = [handles.colCenter-handles.colWidth/2, handles.colCenter+handles.colWidth/2];   
                caxis(handles.axInput, handles.dRange);
                caxis(handles.axOut_1, handles.dRange);
                caxis(handles.axOut_2, handles.dRange);
                caxis(handles.axDiff_1, handles.dRange);
                caxis(handles.axDiff_2, handles.dRange);
        end
    else
        return
    end
catch
end
guidata(hObject, handles);


function CS_LAB_GUI_WindowButtonUpFcn(hObject, eventdata, handles)
% get button release event
handles.FButtonDown = 0;
guidata(hObject, handles)


function pbReset_Callback(hObject, eventdata, handles)
% reset brightness/contrast
handles.dRange = [0,1];
% caxis(handles.axMask, handles.dRange); % just to be safe
caxis(handles.axInput, handles.dRange);
caxis(handles.axOut_1, handles.dRange);
caxis(handles.axOut_2, handles.dRange);
caxis(handles.axDiff_1, handles.dRange);
caxis(handles.axDiff_2, handles.dRange);
guidata(hObject, handles);


function popTrafoDim_2_Callback(hObject, eventdata, handles)
% change sparsifying transformation dimensionality of algo2
iInd = get(hObject,'Value');
if(size(handles.dImg,3) == 1) % just 2D image
    iImgDim = 2;
else % 3D image
    iImgDim = 3;
end
if(get(handles.edSparse,'Value') <= 3)
    iSparseDim = 2; % 2D sparse trafo necessary
else
    iSparseDim = 1; % just 1D sparse trafo necessary
end
if(iInd > iImgDim)
    set(hObject,'Value',iImgDim);
end
if(iInd < iSparseDim)
    set(hObject,'Value',iSparseDim);
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popTrafoDim_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popTrafoDim_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popTrafoDim_1.
function popTrafoDim_1_Callback(hObject, eventdata, handles)
% change sparsifying transformation dimensionality of algo1
iInd = get(hObject,'Value');
if(size(handles.dImg,3) == 1) % just 2D image
    iImgDim = 2;
else % 3D image
    iImgDim = 3;
end
if(get(handles.edSparse,'Value') <= 3)
    iSparseDim = 2; % 2D sparse trafo necessary
else
    iSparseDim = 1; % just 1D sparse trafo necessary
end
if(iInd > iImgDim)
    set(hObject,'Value',iImgDim);
end
if(iInd < iSparseDim)
    set(hObject,'Value',iSparseDim);
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popTrafoDim_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popTrafoDim_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close CS_LAB_GUI.
function CS_LAB_GUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to CS_LAB_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete('parameter_currRun_1.m', 'parameter_currRun_2.m');
delete(hObject);


function popShowMask_Callback(hObject, eventdata, handles)
% switch between mask and PSF
iInd = get(hObject,'Value');

switch iInd
    case 1 % mask
        plotImg = handles.dMask;        
        dRange = [0 1];
    case 2 % psf
        plotImg = handles.dMask;
        fftdim = 1:ndims(plotImg);
        for i=1:length(fftdim)
            plotImg = fftshift(fft(ifftshift(plotImg,fftdim(i)),[],fftdim(i)),fftdim(i));
        end   
        plotImg = abs(plotImg);
        dRange = [0, 0.25*max(plotImg(:))];
end
handles.hPlot.dMask = fPlotAxes(plotImg, 1, handles.hPlot.dMask, handles.axMask, dRange);
guidata(hObject, handles);



function pb_LoadData_Callback(hObject, eventdata, handles)
% load new datasets
[FileName,PathName] = uigetfile('*.mat','Select result file');
if(isempty(FileName))
    return;
end
[~, sFile, sExt] = fileparts(FileName);
handles.cDatasets(:,end+1) = {sFile; [PathName,filesep,FileName]};
iIndDataset = size(handles.cDatasets,2);
set(handles.popDatasets,'String',handles.cDatasets(1,:));
set(handles.popDatasets,'Value',iIndDataset);
[handles.dKSpace, handles.dKSpaceShow, handles.dImg, handles.dMask, handles.measPara, sInfoTxt] = fLoadDataset(handles.cDatasets{2,iIndDataset}, get(handles.edSparse,'Value'));
if(~isempty(handles.dMask)) % mask was supplied
    handles.measPara.prospectiveSub = true;
    handles.dimMask = size(handles.dMask);
    set(handles.edSparse,'Value',handles.measPara.iIndSparseDim);   
    set(handles.edSparse,'Enable','off');
    iStart = 2;
else
    handles.measPara.prospectiveSub = false;
    set(handles.edSparse,'Enable','on');
    iStart = 1;
end
if(size(handles.dImg,3) == 1) % just 2D image
    set(handles.popTrafoDim_1,'Value',2);
    set(handles.popTrafoDim_2,'Value',2);
else % 3D image
    set(handles.popTrafoDim_1,'Value',3);
    set(handles.popTrafoDim_2,'Value',3);
end

% get shown parameter
iInd = get(handles.popShowInput, 'Value');
switch iInd
    case 1 % image
        handles.dInput = scaleImg(handles.dImg, handles.dRange);        
    case 2 % kspace
        handles.dInput = scaleImg(handles.dKSpaceShow, handles.dRange);    
end
handles.iCurrSlice = round(size(handles.dInput,3)/2);

for iI=1:size(handles.idnames,2)
    if(iI > iStart)
        handles.(handles.idnames{2,iI}) = zeros(size(handles.dInput,1),size(handles.dInput,2));
    end
    handles.hPlot.(handles.idnames{2,iI}) = fPlotAxes(handles.(handles.idnames{2,iI}), handles.iCurrSlice, handles.hPlot.(handles.idnames{2,iI}), handles.(handles.idnames{1,iI}), handles.dRange);
end
set(handles.txtInfo,'String',sInfoTxt{1});
set(handles.txtMask,'String',sInfoTxt{2});
guidata(hObject, handles);        
