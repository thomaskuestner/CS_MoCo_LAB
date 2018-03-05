function [image, imageCha, objOut] = CS_reconstruction(path, para)
% CS_RECONSTRUCTION main function to select the appropriate Compressed
% Sensing reconstruction
%
% input: 
% path      - path to file without specifying file extension
%           - path to directory containing files
%           - empty (for file dialog)
%           - kSpace to be reconstructed with dimensionalities:
%             cell array: NSlice x NChannels x NRepetitions x NAverages
%             each cell: 2D: yPhase x xFreq x (nTime/nPhases)
%                        3D: yPhase x xFreq x zPhase
%                        4D: yPhase x xFreq x zPhase x nTime/nPhases
% para      - struct containing to be modified reconstruction parameters
%           and/or struct with handover parameter
%           - path to parameter file
%
%
% output:
% image     reconstructed and sum-of-squared channel combined image
% imageCha  reconstructed individual coil images
% objOut    reconstruction parameters
%
%
% (c) copyright 2012 - 2016 under BSD license
% Thomas Kuestner, University of Tuebingen and University of Stuttgart, Germany
% (thomas.kuestner@{med.uni-tuebingen.de,iss.uni-stuttgart.de})
% 
% Please see license files (./licenses) for accompanied codes (which also 
% include own modifications) of:
% - ESPIRiT v0.1.8  in @ESPIRiT, utils/utils_ESPIRiT     under LICENSE_ESPIRiT       
% - bart v0.2.09    in utils/utils_BART                  under LICENSE_BART
% - SPIRiT v0.3     in @SPIRiT, utils/utils_GRAPPA, 
%                   utils/utils_SPIRiT                   under LICENSE_SPIRiT
% - sparseMRI v0.2  in @sparseMRI, @p2DFT, @A_operator,
%                   @TVOP                                under LICENSE_SPARSEMRI
% - L1Magic v1.11   in @L1_Magic                         under LICENSE_L1MAGIC
% - GIST            in utils/utils_Proximal/GIST         under LICENSE_GIST
% - Wavelab850      in @Wavelet, utils/utils_Wavelet     under LICENSE_WAVELAB
% - Rice Wavelet v3.0 in utils/utils_WaveletRice         under LICENSE_WAVELETRICE
% - Curvelab v2.1.3 in utils/utils_TRAFO/CurveLab-2.1.3  under LICENSE_CURVELAB
% - NUFFT           in @NUFFT, utils/utils_NUFFT         under LICENSE_NUFFT
% - DTCWT           in utils/utils_TRAFO/DTCWT           under LICENSE_DTCWT
% - Shearlet        in utils/utils_TRAFO/Shearlet        under LICENSE_SHEARLET
% - Denoising       in utils/utils_Proximal/Denoising    under LICENSE_DENOISING

warning('off','MATLAB:dispatcher:nameConflict');
warning('off','MATLAB:mat2cell:TrailingUnityVectorArgRemoved');

% global prop

% prepare m-files
currpath = fileparts(mfilename('fullpath'));
addpath(genpath([currpath,filesep,'utils']));
addpath(genpath([currpath,filesep,'io']));
addpath(genpath([currpath,filesep,'preproc']));
addpath(genpath([currpath,filesep,'postproc']));

% load ADCs
if(nargin == 0 || ~exist('path','var') || isempty(path))
    if(strcmp(getenv('computername'),'M123PC'))
        defpath = 'K:\Compressed Sensing\Rohdaten\3D';
    else
        defpath = 'S:\Rohdaten\Thomas\ktFOCUSS\3D';
    end
    [file, path] = uigetfile({'*.mat;*.dat;', 'ADC measdata (*.mat, *.dat)'; '*.h5', 'ISMRMD (*.h5)'} ,'Select measdata file', 'MultiSelect', 'off', defpath);
    if(isequal(file,0))
        error('CS_reconstruction(): User Termination');
    end
    [~,filename,ext] = fileparts(file);
    lReadFromFile = true;
    % in general: [dData, iLC] = fMeasRead(datapath, 'Set', 0, 'Smp', [21 inf]); 
else
    if(ischar(path))
        if(isdir(path))
            [file, path] = uigetfile({'*.mat;*.dat;', 'ADC measdata (*.mat, *.dat)'; '*.h5', 'ISMRMD (*.h5)'} ,'Select measdata file', 'MultiSelect', 'off', path);
            if(isequal(file,0))
                error('CS_reconstruction(): User Termination');
            end
            [~,filename,ext] = fileparts(file);
        else
            [path, filename, ext] = fileparts(path);
        end
        lReadFromFile = true;
    else
        lReadFromFile = false;
    end
end

% read inputs
if(lReadFromFile)
    if(strcmp(ext,'h5')) % ISMRMD
        data = h5read([path,filesep,filename,ext],'dataset/data');
        % convert to processable kSpace
        [ kSpace, measPara ] = convertKSpace_main( data );
    else % ADC measdata (Siemens)
        if(exist(fullfile(path,[filename,'.mat']),'file'))
            matfile = whos('-file', fullfile(path,[filename,'.mat']));
        end
        if(exist(fullfile(path,[filename,'.mat']),'file') && ~ismember('iLC', {matfile.name}))  
            % load mat file to check if it contains iLC/iSP or something else
            data = load(fullfile(path,[filename,'.mat']));
            kSpace = data.kSpace;
            measPara.dim = zeros(1,5);
            measPara.LCall = zeros(1,4);
            dimNames = {'y (phase)', 'x (frequency)', 'z (phase)', 't (time)', 'channels'};
            LCallNames = {'slices', 'acquisitions', 'echos', 'repetitions'};
            dimDefault = ones(1,5);
            LCallDefault = [size(kSpace,1), size(kSpace,4), size(kSpace,5), size(kSpace,3)];
            tmp = size(kSpace{1});
            if(LCallDefault(1) > 1 && length(tmp) > 2)
                dimDefault(1:2) = tmp(1:2);
                dimDefault(4) = tmp(3);
            else
                dimDefault(1:length(tmp)) = tmp;
            end
            dimDefault(5) = size(kSpace,2);
            if(isfield(data,'dim') && isfield(data,'LCall'))
                measPara.dim = data.dim;
                measPara.LCall = data.LCall;
            elseif(isfield(data,'dim') && ~isfield(data,'LCall'))
                measPara.dim = data.dim;
                fprintf('Please specify amount of:\n');
                for i=1:length(measPara.LCall)
                    helper = input(sprintf('%s [%d]: ',LCallNames{i}, LCallDefault(i)),'s');
                    if(isempty(helper)), helper = LCallDefault(i); end;
                    measPara.LCall(i) = helper;
                end

            elseif(~isfield(data,'dim') && isfield(data,'LCall'))
                measPara.LCall = data.LCall;
                fprintf('Please specify dimensions:\n');
                for i=1:length(measPara.dim)
                    helper = input(sprintf('%s [%d]: ',dimNames{i}, dimDefault(i)),'s');
                    if(isempty(helper)), helper = dimDefault(i); end;
                    measPara.dim(i) = helper;
                end

            else
                fprintf('Please specify dimensions:\n');
                for i=1:length(measPara.dim)
                    helper = input(sprintf('%s [%d]: ',dimNames{i}, dimDefault(i)),'s');
                    if(isempty(helper)), helper = dimDefault(i); end;
                    measPara.dim(i) = helper;
                end
                fprintf('Please specify amount of:\n');
                for i=1:length(measPara.LCall)
                    helper = input(sprintf('%s [%d]: ',LCallNames{i}, LCallDefault(i)),'s');
                    if(isempty(helper)), helper = LCallDefault(i); end;
                    measPara.LCall(i) = helper;
                end
            end
            clear 'LCallDefault' 'dimDefault';
            loadedVars = {'dimension', 'oversampling', 'measPara', 'iLC', 'iLCPositions', 'drecksMDH'};
            for i=1:length(loadedVars)
                if(isfield(data,loadedVars{i}))
                    if(any(strcmp(loadedVars{i},{'dimension','oversampling'})))
                        eval(sprintf('measPara.%s = data.%s;', loadedVars{i}, loadedVars{i}));
                    else
                        eval(sprintf('%s = data.%s;', loadedVars{i}, loadedVars{i}));
                    end
                else
                    if(strcmp(loadedVars{i},'measPara'))
                        continue;
                    end
                    eval(sprintf('%s = [];', loadedVars{i}));
                end
            end
            clear 'data';
        else
            [dData, iLC, iEvalInfoMask, drecksMDH, iLCPositions] = fMeas_main(fullfile(path, [filename, ext]));
            measPara = [];
        end
    end
else
    % kSpace was directly handed over
    kSpace = path;
    path = currpath; % for compatibility later
    filename = '';
    if(isstruct(para))
        loadedVars = {'measPara', 'dimension', 'oversampling', 'iLC', 'iLCPositions', 'drecksMDH', 'savePath'};
        for i=1:length(loadedVars)
            if(isfield(para,loadedVars{i}))
                if(any(strcmp(loadedVars{i},{'dimension','oversampling'})))
                    eval(sprintf('measPara.%s = para.%s;', loadedVars{i}, loadedVars{i})); 
                else
                    eval(sprintf('%s = para.%s;', loadedVars{i}, loadedVars{i}));
                end
                eval(sprintf('para = rmfield(para,''%s'');', loadedVars{i}));
            else
                eval(sprintf('%s = [];', loadedVars{i}));
            end
        end
    end
end
    
%% save output
prop.flagSave = false;


%% show output in imagine
prop.flagPlot = true;


%% display output
prop.disp.names = {'Repetitions', 'Averages', 'Slices', 'Channels', 'Time', 'Frequency', 'Extracting K-Space', 'Phase Correction', 'GRAPPA calibration', 'FOCUSS', 'CG', 'ESPReSSo', 'Trafo', 'POCS', 'MC', 'Line Search', 'Log barrier', 'Newton', 'Proximal Average', 'ADMM', 'SENSE - Slice', 'SENSE - Time', 'SENSE - Calibration', 'SENSE - Kernel', 'SENSE - Concatenate'};
if(usejava('jvm') && ~feature('ShowFigureWindows'))
    % use text-based alternative
    prop.flagDisp = false;
    prop.flagPlot = prop.flagPlot & false;
    
    % initialize some console variables
    prop.maxvals = zeros(1,length(prop.disp.names));
    prop.lastOut = '';
else
    % use GUI dialogs
    prop.flagDisp = true;
    prop.flagPlot = prop.flagPlot & true;
    
    % initialize some GUI variables
    prop.disp.openBar = false(1,size(prop.disp.names,2));
    prop.disp.colors = colorGradient(rgb('DarkRed'),rgb('ForestGreen'),100);
    prop.disp.ppm = cell(1,length(prop.disp.names));
    % initialize some console variables
    prop.maxvals = zeros(1,length(prop.disp.names));
    prop.lastOut = '';
end


%% parallel computing
prop.flagParallel = false;
prop.openPool = false;
helper = ver;
lFound = false;
for i=1:length(helper)
    if(strcmp(helper(i).Name,'Parallel Computing Toolbox'))
        lFound = true;
        prop.flagParallel = prop.flagParallel & true;
        if(verLessThan('matlab','8.4'))
            if(matlabpool('size') > 0)
    %             matlabpool close
                prop.openPool = true;
            end
        else
            if(~isempty(gcp('nocreate')))
                prop.openPool = true;
            end
        end
        break;
    end
end
if(~lFound), prop.flagParallel = false; end
    


%% load reconstruction parameter
if(nargin > 1 && exist('para', 'var') && ischar(para))
    [pathPara, filenamePara, ~] = fileparts(para);
    cd(pathPara);
    eval([filenamePara,';']);
    cd(currpath);
else
    % load default parameter set
    eval('parameters_default;'); 
end


%% replace parameters (if necessary)
if(nargin > 1 && exist('para', 'var') && isstruct(para))
    names = fieldnames(para);
    for i=1:length(names)
        if(isstruct(para.(names{i})))
            innerNames = fieldnames(para.(names{i}));
            for j=1:length(innerNames)
                eval(sprintf('%s.(innerNames{j}) = para.(names{i}).(innerNames{j});',names{i}));
                if(strcmp(innerNames{j},'trafodim') && ~any(strcmp(innerNames,'scrambledim')))
                    eval(sprintf('%s.scrambledim = para.(names{i}).(innerNames{j})(1,:);',names{i}));
                end
            end
            clear 'innerNames'
        else
            eval(sprintf('%s = para.%s;',names{i},names{i}));
        end
    end
    clear 'names'
end
clear 'para';
% save MC parameters
if(exist('savePath','var'))
    mc.savePath = savePath;
else
    mc.savePath = '';
end
% check for missing necessary parameters
loadedVars = {'measPara', 'iLC', 'iLCPositions', 'drecksMDH'};
for iI=1:length(loadedVars)
    if(eval(sprintf('~exist(''%s'',''var'')',loadedVars{iI})))
        eval(sprintf('%s = [];', loadedVars{iI}));
    end
end


%% initialize dispProgress
dispProgress('InitDispProgress',prop);


%% get data dimension and loop counters
[measPara, espresso, postproc] = extractMeasPara( iLC, measPara, espresso, postproc, flagOversampling, drecksMDH, iLCPositions );


%% check variable consistency and existence
allPara = { cstype,   measPara.dimension,   measPara.dim,   iNOUTER,   iNINNER,   p,   lambda,   epsilon,   measPara.oversampling,   lambdaCalib,   calibTyk,   reconTyk,   kernelSize,   calibSize,   opt_problem,   postproc,   trafo,   FFTwindow,   espresso,   measPara.LCall, lbtol,   tvmu,   cgtol,   cgmaxiter,   flagToolbox,   flagZeropadding,  measPara.aniso,  measPara.precision, reconDIM,   flags,   lambdaGroup,   lambdaNLTV,   lambdaTV;...
           'cstype', 'dimension',           'dim',         'iNOUTER', 'iNINNER', 'p', 'lambda', 'epsilon', 'oversampling',          'lambdaCalib', 'calibTyk', 'reconTyk', 'kernelSize', 'calibSize', 'opt_problem', 'postproc', 'trafo', 'FFTwindow', 'espresso', 'LCall',        'lbtol', 'tvmu', 'cgtol', 'cgmaxiter', 'flagToolbox', 'flagZeropadding', 'aniso',        'precision',        'reconDIM', 'flags', 'lambdaGroup', 'lambdaNLTV', 'lambdaTV'}; %#ok<*NODEF>

[consExist, msg, msgWarn, kernelSize, calibSize, FFTwindow, trafo, flagToolbox, measPara.precision, flags, reconDIM] = checkConsExist(allPara);
if(~consExist)
    error('CS_reconstruction:checkConsExist','%s', msg); 
end
if(~isempty(msgWarn))
    warning('CS_reconstruction:checkConsExist\n%s', msgWarn);
end
clear 'allPara' 'msg' 'msgWarn' 'consExist';
            

%% initialize recon objects
if(strcmp(cstype,'FOCUSS'))
    obj = FOCUSS(iNOUTER, iNINNER, p, lambda, epsilon, measPara, lambdaCalib, lambdaTV, lambdaESPReSSo, lambdaMC, kernelSize, calibTyk, calibSize, flagZeropadding, FFTwindow, trafo, espresso, sScaling, mc, lSense);
        
elseif(strcmp(cstype,'SPIRiT_CG') || strcmp(cstype,'SPIRiT_POCS'))
    obj = SPIRiT_Wrapper(cstype(8:length(cstype)), iNINNER, measPara, lambda, kernelSize, calibTyk, reconTyk, trafo, espresso, FFTwindow, calibSize); %#ok<*COLND>
        
elseif(strcmp(cstype,'ESPIRiT_CG') || strcmp(cstype,'ESPIRiT_L1'))
    obj = ESPIRiT_Wrapper(cstype(9:length(cstype)), iNINNER, iNIterSplit, measPara, kernelSize, calibTyk, reconTyk, eigThresh_k, eigThresh_im, n_maps, splitWeight, trafo, flagToolbox, path, espresso, FFTwindow, calibSize);

elseif(strcmp(cstype,'BART'))
    obj = BART_Wrapper(lambda, lambdaCalib, lambdaTV, lambdaMC, measPara, kernelSize, n_maps, trafo, espresso, FFTwindow, currpath, mc, calibSize);
    
elseif(strcmp(cstype,'sparseMRI'))
    obj = sparseMRI(iNOUTER, iNINNER, measPara, lambda, lambdaTV, lambdaESPReSSo, p, trafo, espresso, FFTwindow, l1Smooth, lineSearchItnlim, lineSearchAlpha, lineSearchBeta, lineSearchT0, gradToll);
    
elseif(strcmp(cstype,'L1_Magic_TV') || strcmp(cstype,'L1_Magic_L1') || strcmp(cstype,'L1_Magic_TVDantzig') || strcmp(cstype,'L1_Magic_L1Dantzig'))
    obj = L1_Magic(cstype(10:length(cstype)), measPara, epsilon, lbtol, tvmu, cgtol, cgmaxiter, pdmaxiter, espresso);

elseif(strcmp(cstype,'Zero'))
    obj = Zero(measPara, espresso);
    
elseif(strcmp(cstype,'FCSA_WaTMRI') || strcmp(cstype,'FCSA_SLEP') || strcmp(cstype,'FCSA_proxA') || strcmp(cstype,'BFCSA_proxA') || strcmp(cstype,'ADMM_proxA') || strcmp(cstype,'SB_proxA') || ...
        strcmp(cstype,'A_PFISTA') || strcmp(cstype,'A_TDIHT') || strcmp(cstype,'A_ADMM_L1') || strcmp(cstype,'A_ADMM_L0') || strcmp(cstype,'A_ADMM_SCAD') || strcmp(cstype,'A_ADMM_MCP') || strcmp(cstype,'A_ADMM_ATAN') || strcmp(cstype,'A_ADMM_PSHRINK'))
    obj = Proximal(cstype, iNINNER, itrNLTV, iNOUTER, measPara, trafo, lambda, lambdaTV, lambdaGroup, lambdaNLTV, lambdaNLTV_h, mue, regularizerWeights, flags, reconDIM, espresso);
    
else
    error('CS_reconstruction(): Undefined reconstruction cstype');
end


%% optimize fft command
if(prop.flagOptim)
    optimizeFFT(measPara.dimension, measPara.dim, kernelSize, fft_planner_method, trafo.fftdim, lambdaCalib, flagZeropadding, measPara.precision);
end


%% extract kSpace data from measurement file
if(exist('iLC','var') && ~isempty(iLC) && ~exist('kSpace', 'var'))
    [obj.kSpace, evalMask] = extractKSpace_main(dData, iLC, iEvalInfoMask, measPara);
else
    obj.kSpace = kSpace;
    clear 'kSpace';
end


%% correct phase perturbations (echo alignment)
if(isfield(measPara,'sequenceName') && strcmp(measPara.sequenceName, 'CS_EPI'))
    obj.kSpace = alignEchos(obj.kSpace, evalMask, measPara, dData, iLC, iEvalInfoMask);
end


% clear all unnecessary variables
clear 'dData' 'iLC' 'dDataRep' 'iLCRep' 'dDataSli' 'iLCSli' 'drecksMDH' 'iLCPositions' 'evalMask' 'helper' 'espresso' 'trafo' 'flagOversampling' 'mc';
sCallingStack = dbstack;
if(strcmp(sCallingStack(end).name,'fRetroGateAndRecon') || strcmp(sCallingStack(end).name,'MotionCorrection')) % for 4D CS_Retro
    evalin('caller', 'clear ''kSpace'';'); % ''para''
end
clear 'sCallingStack';


%% start reconstruction
cs_time = tic;
imageCha = cell(size(obj.kSpace));

if(lAutoCheck && ( ((strcmp(measPara.dimension,'3D') || strcmp(measPara.dimension,'4D')) && nnz(obj.kSpace{1,1}) >= 0.95*prod(measPara.dim(1:4))) || ...
        (strcmp(measPara.dimension,'2D') && nnz(obj.kSpace{1,1}) >= 0.95*prod(measPara.dim([1:2,4]))) || ...
        (obj.espresso.state && ~isempty(regexp(filename,'\w*x1\D\w*','once'))) || (obj.espresso.state && isfield(measPara, 'CSAcceleration') && measPara.CSAcceleration == 1)))
    % check for acceleration factor == 1 ==> full acquisition
    if(strcmp(measPara.dimension,'2D'))
        fftdim = 1:2;
        pfdim = [3 1 2];
        subs = {':', ':', ':', ':'}; % y-x-cha
        subsid = 4;
        if(measPara.dim(4) > 1)
            pfdim = [4 1 2 3];
            subs = {':', ':', ':', ':'}; % y-x-t-cha
            subsid = 3;
        end
    elseif(strcmp(measPara.dimension,'3D'))
        fftdim = 1:3;
        pfdim = [4 1 2 3];
        subs = {':', ':', ':', ':', ':'}; % y-x-z-cha
        subsid = 5;
    elseif(strcmp(measPara.dimension,'4D'))
        fftdim = 1:3;
        pfdim = [5 1 2 3 4];
        subs = {':', ':', ':', ':', ':'}; % y-x-z-t-cha
        subsid = 4;
    end
    totalCount = measPara.LCall(2)*measPara.LCall(4)*measPara.LCall(1)*measPara.dim(4);
    dispProgress('Repetitions', 0, measPara.LCall(4));
    dispProgress('Averages', 0, totalCount);
    dispProgress('Time', 0, totalCount);
    for iRep=1:measPara.LCall(4)
        for iAvg=1:measPara.LCall(2)
            if(obj.espresso.state)
                for iSli=1:measPara.LCall(1)
                    kSpaceTmp = squeeze(cell2mat(shiftdim(obj.kSpace(iSli,:,iRep,iAvg),-2)));
                    imgTmp = zeros(size(kSpaceTmp));                 
                    for iTime=1:measPara.dim(4)
                        subs{subsid} = iTime;
                        imgTmp(subs{:}) = ipermute(espresso_recon(permute(kSpaceTmp(subs{:}),pfdim), obj.espresso.iter, false),pfdim);
                        dispProgress('Time', ((iRep-1)*measPara.LCall(2)*measPara.LCall(1)*measPara.dim(4) + (iAvg-1)*measPara.LCall(1)*measPara.dim(4) + (iSli-1)*measPara.dim(4) +iTime)/totalCount);
                    end
                    if(strcmp(measPara.dimension,'2D') && measPara.dim(4) == 1) % 2D
                        imageCha(iSli,:,iRep,iAvg) = shiftdim(squeeze(mat2cell(imgTmp, size(obj.kSpace{1,1},1),size(obj.kSpace{1,1},2),ones(1,measPara.dim(5)))),-1);
                    else % 2Dt, 3D, 4D
                        imageCha(iSli,:,iRep,iAvg) = shiftdim(squeeze(mat2cell(imgTmp, size(obj.kSpace{1,1},1),size(obj.kSpace{1,1},2),size(obj.kSpace{1,1},3),ones(1,measPara.dim(5)))),-1);
                    end
                    clear 'imgTmp' 'kSpaceTmp';
                end
            else              
                imageCha(:,:,iRep,iAvg) = cellfun(@(x) ifftnshift(x,fftdim), obj.kSpace(:,:,iRep,iAvg), 'UniformOutput', false);
            end  
            dispProgress('Averages', ((iRep-1)*measPara.LCall(2)*measPara.LCall(1)*measPara.dim(4) + iAvg*measPara.LCall(1)*measPara.dim(4))/totalCount);
        end
        dispProgress('Repetitions', iRep/measPara.LCall(4));
    end
    dispProgress('Time','Close');
else
    % Compressed Sensing reconstruction
    dispProgress('Repetitions', 0, measPara.LCall(4));
    dispProgress('Averages', 0, measPara.LCall(4)*measPara.LCall(2));
    for iRep=1:measPara.LCall(4)
        for iAvg=1:measPara.LCall(2)
            imageCha(:,:,iRep,iAvg) = obj.recon_main(iRep,iAvg); % obj.kSpace(:,:,iRep,iAvg)
            dispProgress('Averages', ((iRep-1)*measPara.LCall(2) + iAvg)/(measPara.LCall(4)*measPara.LCall(2)));
        end
        dispProgress('Repetitions', iRep/measPara.LCall(4));
    end
end 
if(exist('cutOutMask','var')), imageCha = imageCha(cutOutMask); end;
dispProgress('Averages', 'Close');
dispProgress('Repetitions', 'Close');
obj.kSpace = []; % clear kSpace variable

cs_time = toc(cs_time);
dRecontime(1) = mod(cs_time,60);
dRecontime(2) = mod((cs_time - dRecontime(1))/60,60);
dRecontime(3) = mod(((cs_time - dRecontime(1))/60 - dRecontime(2))/60, 24);
if(cs_time >= 3600)   
    fprintf('Execution time: %.2dh %.2dmin %.2ds\n', dRecontime(3), dRecontime(2), round(dRecontime(1)));
elseif(cs_time >= 60)
    fprintf('Execution time: %.2dmin %.2ds\n', dRecontime(2), round(dRecontime(1)));
else
    fprintf('Execution time: %.2ds\n', round(dRecontime(1)));
end
obj.dRecontime = dRecontime;

%% postprocessing
if(isempty(obj.kernelSize))
    obj.kernelSize = kernelSize;
end

image = cell(1,measPara.LCall(4));
for iRep=1:measPara.LCall(4)
    imageAvg = cell(1,measPara.LCall(2));
    for iAvg=1:measPara.LCall(2)
        % each element: 2D (y-x + one slice), 3D (y-x-t + one slice, y-x + several slices, y-x-z) or
        % 4D (y-x-t + several slices, y-x-z-t)
        imageAvg{1,iAvg} = postproc_main(obj, imageCha(:,:,iRep,iAvg), postproc);
    end
    image{iRep} = squeeze(mean(cell2mat(shiftdim(imageAvg,-4)),6));
end
image = squeeze(cell2mat(shiftdim(image,-3)));

% % average over repetitions
% image = squeeze(mean(cell2mat(shiftdim(image,-3)),5));


%% delete fftw_wisdom variable
if(prop.flagOptim)
    evalin('base', 'clear ''fftw_wisdom''');
end

%% save image
if(nargout > 2 || prop.flagSave)
    % convert obj->struct and reduce size
    objProp = properties(obj);
    for iObj = 1:numel(objProp)
        if(any(strcmp(objProp{iObj},{'kernel','kernelImg','kCalib','A','kSpace'})))
            continue;
        end
        eval(['objOut.',objProp{iObj},' = obj.',objProp{iObj},';']);
    end
end
if(prop.flagSave)
    if(exist('filename','var') && ~isempty(filename))
        savename = [filename,'_',cstype,'_',transformation];
    else
        savename = [cstype,'_',transformation];
    end
    dirContent = dir(path);
    occurence = 0;
    for i=1:length(dirContent)
        if(dirContent(i).isdir)
            continue;
        end
        tmp = regexp(dirContent(i).name,[savename,'[_]*\d*'],'once');
%         if(strcmpi(dirContent(i).name,savename))
        if(~isempty(tmp))
            occurence = occurence + 1;
        end
    end
    if(occurence > 0)
         savename = [savename,'_',num2str(occurence)];
    end
%     save([path,filesep,savename,'.mat'], 'image', 'imageCha', 'objOut');
    save([path,filesep,savename,'.mat'], 'image', 'objOut', '-v7.3');
end


%% plot resulting image
if(prop.flagPlot && exist('imagine','file'))
    if(iscell(image))
        imagePlot = cell2mat(shiftdim(image,-1));
    else
        imagePlot = image;
    end
    
    % set displayed name
    figname = [cstype,'_',transformation];
    tmp = functions(FFTwindow.type);
    if(~strcmp(tmp.function,'rectwin'))
        figname = [figname,'_',tmp.function];
    end    
    if(obj.espresso.state)
        figname = [figname,'_espresso'];
    end
    figname = [figname,'_',postproc.type];
    
    % if imagine is already open, add new image in same window
    hfig = findobj('type','figure','name','IMAGINE 2.0 Belly Jeans');
    if(length(hfig) > 1)
        hfig = max(hfig);
    end
    
    if(~isempty(hfig) && ishandle(hfig))
        fh = get(hfig,'userdata');
        data = fh();
        
        
        oldImage = cell(1,size(data,2));
        oldTitle = oldImage;
        cmd = '';
        for i=1:size(data,2)
            oldImage{i} = data(i).dImg;
            oldTitle{i} = data(i).sName;
            cmd = [cmd, sprintf('oldImage{%d},''Name'',oldTitle{%d},',i,i)];
        end

        close(hfig);
        cmd = ['imagine(',cmd, 'imagePlot,''Name'',figname);'];
        eval(cmd);
        clear oldImage oldTitle cmd hfig;
    else
        imagine(imagePlot,'Name',figname);
    end
end


end