function [consExist, msg, msgWarn, kernelSize, calibSize, window, trafo, flagToolbox, sPrecision, flags, reconDIM] = checkConsExist(allPara)
% check consistency and existence for all needed variables
% allPara = {type, dimension, dim, iNOUTER, iNINNER, 
%            p, lambda, epsilon, oversampling, lambdaCalib, 
%            calibTyk, reconTyk, kernelSize, calibSize, opt_problem, 
%            postproc, trafo, FFTwindow, espresso, LCall, 
%            lbtol, tvmu, cgtol, cgmaxiter, flagToolbox,
%            flagToolbox, aniso, precision, reconDim, flags, 
%            lambdaGroup, lambdaNLTV, lambdaTV}
%
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

consExist = true;

allVarnames = allPara(2,:);
allPara = allPara(1,:);
% get positions of variables in allPara cell array
pcstype = find(strcmp('cstype',allVarnames));
pdimension = find(strcmp('dimension',allVarnames));
pdim = find(strcmp('dim',allVarnames));
piNOUTER = find(strcmp('iNOUTER',allVarnames));
piNINNER = find(strcmp('iNINNER',allVarnames));
pp = find(strcmp('p',allVarnames));
plambda = find(strcmp('lambda',allVarnames));
pepsilon = find(strcmp('epsilon',allVarnames));
poversampling = find(strcmp('oversampling',allVarnames));
plambdaCalib = find(strcmp('lambdaCalib',allVarnames));
pcalibTyk = find(strcmp('calibTyk',allVarnames));
preconTyk = find(strcmp('reconTyk',allVarnames));
pkernelSize = find(strcmp('kernelSize',allVarnames));
pcalibSize = find(strcmp('calibSize',allVarnames));
popt_problem = find(strcmp('opt_problem',allVarnames));
ppostproc = find(strcmp('postproc',allVarnames));
ptrafo = find(strcmp('trafo',allVarnames));
pwindow = find(strcmp('FFTwindow',allVarnames));
pespresso = find(strcmp('espresso',allVarnames));
plbtol = find(strcmp('lbtol',allVarnames));
pmu = find(strcmp('tvmu',allVarnames));
pcgtol = find(strcmp('cgtol',allVarnames));
pcgmaxiter = find(strcmp('cgmaxiter',allVarnames));
pflagToolbox = find(strcmp('flagToolbox',allVarnames));
pflagZeropadding = find(strcmp('flagZeropadding',allVarnames));
paniso = strcmp('aniso',allVarnames);
pPrecision = find(strcmp('precision',allVarnames));
preconDim = find(strcmp('reconDIM',allVarnames));
pflags = find(strcmp('flags',allVarnames));
pLambdaGroup = find(strcmp('lambdaGroup',allVarnames));
pLambdaNLTV = find(strcmp('lambdaNLTV',allVarnames));
pLambdaTV = find(strcmp('lambdaTV',allVarnames));

kernelSize = allPara{1,pkernelSize};
calibSize = allPara{1,pcalibSize};
cstype = allPara{1,pcstype};
dimension = allPara{1,pdimension};
dim = allPara{1,pdim};
aniso = allPara{1,paniso};
window = allPara{1,pwindow};
trafoType = allPara{1,ptrafo}.trafoType;
trafo = allPara{1,ptrafo};
flagToolbox = allPara{1,pflagToolbox};
sPrecision = allPara{1,pPrecision};
flags = allPara{1,pflags};
reconDIM = allPara{1,preconDim};

emptyPara = cellfun(@isempty, allPara(1,:));
range = [];
necTrafoPara = [];
msg = '';
msgWarn = '';

% check type, dimension and postproc
if(~strcmp(cstype,'FOCUSS') && ~strcmp(cstype,'sparseMRI') && ...
        ~strcmp(cstype,'SPIRiT_CG') && ~strcmp(cstype,'SPIRiT_POCS') && ...
        ~strcmp(cstype,'ESPIRiT_CG') && ~strcmp(cstype,'ESPIRiT_L1') && ...
        ~strcmp(cstype,'L1_Magic_TV') && ~strcmp(cstype,'L1_Magic_L1') && ...
        ~strcmp(cstype,'L1_Magic_L1Dantzig') && ~strcmp(cstype,'L1_Magic_TVDantzig') && ...
        ~strcmp(cstype,'RecPF') && ~strcmp(cstype,'BART') && ...
        ~strcmp(cstype,'FCSA_TV') && ~strcmp(cstype,'FCSA_NLTV') && ...
        ~strcmp(cstype,'FCSA_WaTMRI') && ~strcmp(cstype,'FCSA_SLEP') && ...
        ~strcmp(cstype,'FCSA_proxA') && ~strcmp(cstype,'BFCSA_proxA') && ...
        ~strcmp(cstype,'ADMM_proxA') && ~strcmp(cstype,'SB_proxA') && ...
        ~strcmp(cstype,'A_PFISTA') && ~strcmp(cstype,'A_TDIHT') && ...
        ~strcmp(cstype,'A_ADMM_L1') && ~strcmp(cstype,'A_ADMM_L0') && ...
        ~strcmp(cstype,'A_ADMM_SCAD') && ~strcmp(cstype,'A_ADMM_MCP') && ...
        ~strcmp(cstype,'A_ADMM_ATAN') && ~strcmp(cstype,'A_ADMM_PSHRINK') && ...
        ~strcmp(cstype,'Zero'))
    msg = 'Undefined reconstruction type';
    consExist = false;
    return
end

if(~strcmp(trafoType,'fft') && ~strcmp(trafoType,'wavelet_mat') && ...
        ~strcmp(trafoType,'wavelet_lab') && ~strcmp(trafoType,'mellin') && ...
        ~strcmp(trafoType,'pca') && ~strcmp(trafoType,'dct') && ...
        ~strcmp(trafoType,'surfacelet') && ~strcmp(trafoType,'curvelet'))
    msg = 'Undefined sparsifying transformation';
    consExist = false;
    return
end

if(~strcmp(dimension,'2D') && ~strcmp(dimension,'3D') && ~strcmp(dimension,'4D') && ~strcmp(dimension,'5D'))
    msg = 'Undetermined dimensionality';
    consExist = false;
    return
end

if(~strcmp(allPara{1,ppostproc}.type,'rss') && ~strcmp(allPara{1,ppostproc}.type,'SENSE'))
    msg = 'Undefined post reconstruction type';
    consExist = false;
    return
end


%% existence
if(any(emptyPara([pcstype,pdimension,piNOUTER,piNINNER,pp,plambda,pepsilon,ppostproc])))
    % necessary: type, dimension, iNOUTER, iNINNER, p, lambda, epsilon,
    %            postproc
    out = '';
    for i=find(emptyPara([1,2,5:9,16]))
        out = [out, allVarnames{1,i}, ' '];
    end
    msg = sprintf('Not all necessary algorithm parameters are set\nMissing parameter(s): %s', out);
    consExist = false;
    return
end

if(strcmp(cstype,'FOCUSS'))
    % necessary: lambdaCalib, calibTyk, kernelSize
    range = [plambdaCalib,pcalibTyk,pkernelSize];
elseif(strcmp(cstype,'SPIRiT_CG') || strcmp(cstype,'SPIRiT_POCS') || strcmp(cstype,'ESPIRiT_CG') || strcmp(cstype,'ESPIRiT_L1'))
    % necessary: wavWeight, calibTyk, reconTyk, kernelSize
    range = [ptrafo,pcalibTyk,preconTyk,pkernelSize];
elseif(strcmp(cstype,'L1_Magic_TV'))
    % necessary: epsilon, lbtol, mu, cgtol, cgmaxiter
    range = [pepsilon, plbtol, pmu, pcgtol, pcgmaxiter];
end

if(~isempty(range))
    if(any(emptyPara(range)))
        out = '';
        for i=find(emptyPara(range))
            out = [out, allVarnames{1,i}, ' '];
        end
        msg = sprintf('Not all necessary algorithm parameters are set\nMissing parameter(s): %s', out);
        consExist = false;
        return
    end
end

if(strcmp(trafoType,'fft'))
    necTrafoPara = {'windowing', 'windowType', 'windowOpt'};
elseif(strcmp(trafoType,'wavelet_mat'))
    necTrafoPara = {'waveletFilter', 'waveletFilterSize', 'extMode', 'waveletStages', 'flagThresholding', 'wavWeight'};
elseif(strcmp(trafoType,'wavelet_lab'))
    necTrafoPara = {'waveletFilter', 'waveletFilterSize', 'waveletStages', 'wavWeight'};
elseif(strcmp(trafoType,'mellin'))
    necTrafoPara = {'extrapVal', 'trafoSize', 'bTrafoType', 'interp', 'sigma'};
elseif(strcmp(trafoType,'curvelet'))
    necTrafoPara = { 'nbscales', 'allCurvelets', 'nbdstz_coarse'};
elseif(strcmp(trafoType,'surfacelet'))
    necTrafoPara = {'Pyr_mode', 'HGfname', 'bo', 'msize', 'beta', 'lambda', 'decompLevels', 'Lev_array', 'padOnSameSize'};
end

if(~isempty(necTrafoPara))
    out = '';
    for i=1:length(necTrafoPara)
        if(~isfield(allPara{1,ptrafo},necTrafoPara{i}))
            out = [out, necTrafoPara{i}, ' '];
        end
    end
    if(~isempty(out))
        msg = sprintf('Not all necessary algorithm parameters are set\nMissing parameter(s): %s', out);
        consExist = false;
        return
    end
end

if(~emptyPara(pwindow))
    fname = functions(window.type);
    fname = fname.function;
    switch fname
        case {'rectwin', 'barthannwin', 'bartlett', 'bohmanwin', 'parzenwin', 'triang'}
            % nop
        case {'blackman','blackmanharris', 'hamming', 'hann', 'flattopwin', 'nuttallwin', 'lanczos', 'welch', 'cosine'}
            if(isempty(window.windowOpt{1}))
                msg = sprintf('Window option symmetric or periodic not set for cut-out window');
                consExist = false;
                return
            end
        case {'chebwin', 'tukeywin', 'gausswin', 'kaiser', 'planckTaper'}
            if(isempty(window.windowOpt{1}) || isempty(window.windowOpt{2}))
                msg = sprintf('Window option symmetric or periodic and additional window parameter must be set');
                consExist = false;
                return
            end
        otherwise
            msg = sprintf('Unknown cut-out window type');
            consExist = false;
            return
    end
end    

if(~emptyPara(pespresso))
    if(allPara{1,pespresso}.state)
        if(isempty(allPara{1,pespresso}.iter))
            msg = sprintf('Partial Fourier reconstruction chosen, but amount of iterations is not set');
            consExist = false;
        end
    end
end


%% consistency
for i=[piNOUTER,piNINNER,pp,plambda,pepsilon,poversampling,plambdaCalib,pcalibTyk,preconTyk,plbtol,pmu,pcgtol,pcgmaxiter]
    % iNOUTER, iNINNER, p, lambda, epsilon, oversampling, lambdaCalib,
    % calibTyk, lbtol, mu, cgtol, cgmaxiter
    switch i
        case poversampling
            if(~isvector(allPara{1,i}{1,1}) || any(allPara{1,i}{1,1} <= 0) || allPara{1,i}{1,1}(end) > dim(2))
                msg = 'antiAliasing in x-direction must be a vector with elements larger than 0';
                consExist = false;
                return
            end
            if(~isempty(allPara{1,i}{1,2}) && (~isvector(allPara{1,i}{1,2}) || any(allPara{1,i}{1,2} <= 0) || (~isempty(aniso) && allPara{1,i}{1,2}(end) > aniso(1))))
                msg = 'phase oversampling correction must be a vector with elements larger than 0';
                consExist = false;
                return
            end
            if(~isempty(allPara{1,i}{1,3}) && (~isvector(allPara{1,i}{1,3}) || any(allPara{1,i}{1,3} <= 0) || (strcmp(dimension,'2D') && allPara{1,i}{1,3}(end) > allPara{1,20}(1)) || (strcmp(dimension,'3D') && ~isempty(aniso) && allPara{1,i}{1,3}(end) > aniso(3))))
                msg = 'slice oversampling correction must be a vector with elements larger than 0';
                consExist = false;
                return
            end
        otherwise
            if(~isscalar(allPara{1,i}) || allPara{1,i} < 0)
                msg = sprintf('%s must be a positive scalar', allVarnames{1,i});
                consExist = false;
                return
            end
    end
end

% if(~isempty(iTotalLines))
%     if(iTotalLines ~= dim(1))
%         msg = 'set iTotalLines parameter does not fit with calculated iTotalLines parameter';
%         consExist = false;
%         return
%     else
%         iTotalLines = dim(1);    
%     end
% else
%     iTotalLines = dim(1);
% end

% adapt kernelSize and calibSize
if(strcmp(dimension,'2D'))
    if(size(kernelSize,2) > 2)
        kernelSize = [kernelSize(1) kernelSize(3)];
        msgWarn = [msgWarn,'Changing kernelSize to 2D mode; '];
    end
    if(any(kernelSize > dim(1:2)) &&  allPara{1,plambdaCalib} > 0)
        msg = 'kernelSize must be smaller than image dimensions';
        consExist = false;
        return
    end
    if(~isempty(calibSize))
        if(size(calibSize,2) > 2)
            calibSize = [calibSize(1) calibSize(3)];
            msgWarn = [msgWarn,'Changing calibSize to 2D mode; '];
        end
        if(any(calibSize > dim(1:2)) &&  allPara{1,plambdaCalib} > 0)
            msg = 'calibSize must be smaller than image dimensions';
            consExist = false;
            return
        end
    end
    
    % for fftdim and scrambledim check
    if(dim(4) > 1) % time direction t-y-x
        testdim = 3;
        permDim = [4 1 2 3];
    else % y-x
        testdim = 2;
        permDim = [1 2 3 4];
    end        
        
elseif(strcmp(dimension,'3D') || strcmp(dimension,'4D') || strcmp(dimension,'5D'))
    if(size(kernelSize,2) < 3)
        msg = '3-dimensional kernel size needed for 3D/4D patient data sets';
        consExist = false;
        return
    end
%     if(any(kernelSize > [dim(1) dim(3) dim(2)]) &&  (allPara{1,plambdaCalib} > 0 || strcmp(cstype,'SPIRiT_CG') || strcmp(cstype,'SPIRiT_POCS') || strcmp(cstype,'ESPIRiT_CG') || strcmp(cstype,'ESPIRiT_L1')))
%         msg = 'kernelSize must be smaller than image dimensions';
%         consExist = false;
%         return
%     end
    if(~isempty(calibSize))
        if(size(calibSize,2) < 3)
            msg = '3-dimensional calibration size needed for 3D/4D patient data sets';
            consExist = false;
            return
        end
        if(any(calibSize > [dim(1) dim(3) dim(2)]) &&  allPara{1,plambdaCalib} > 0)
            msg = 'calibSize must be smaller than image dimensions';
            consExist = false;
            return
        end
    end
    
    % for fftdim and scrambledim check
    if(strcmp(dimension,'3D')) % y-z-x
        testdim = 3;
        permDim = [1 3 2 4];
    else % t-y-z-x
        testdim = 4;
        permDim = [4 1 3 2];
    end
end

% translate dimensions to be transformed
if(trafo.scrambledim(4))
    msgWarn = [msgWarn,'Time dimension should not be scrambled!; '];
    trafo.scrambledim(4) = 0; % time dimension should not be scrambled
end
% consistency check
if(any(any(trafo.trafodim(:,[1 3 4]) ~= 1 & trafo.trafodim(:,[1 3 4]) ~= 0)) && any(trafo.trafodim(:,2) ~= 0) && any(trafo.trafodim(:,2)) ~= 1 && trafo.trafodim(1,2) ~= 2)
    msg = 'Invalid input of trafo dimensions';
    consExist = false;
    return    
end
% perform fourier transformation before or after reconstruction
% allowing the transformation just to take place between A_cart and A_newBase 
trafo.fftBA = false(2,4);
trafo.fftBA(1,trafo.trafodim(1,permDim) == 2) = true; % 1st line: before CS recon
trafo.fftBA(2,trafo.trafodim(1,permDim) == 0) = true; % 2nd line: after CS recon (except for Time)
trafo.fftBA(2,permDim == 4) = false;
lIdx = trafo.trafodim(1,:) == 2;
trafo.trafodim(1,lIdx) = 0; % to avoid fourier transformation in compBasis (also if kspaceTrafo is chosen)
trafo.scrambledim(1,lIdx) = 0;
if(any(any(trafo.fftBA(:,testdim+1:end))))
    trafo.fftBA(:,testdim+1:end) = false;
end
% % if(trafo.trafodim(1,2) == 2)
% %     trafo.fftX = [true, false]; % [before, after] CS recon
% %     trafo.trafodim(1,2) = 0; 
% % elseif(trafo.trafodim(1,2) == 1)
% %     trafo.fftX = [false, false];
% % elseif(trafo.trafodim(1,2) == 0)
% %     trafo.fftX = [false, true];
% % end
trafo.fftdim = find(trafo.trafodim(1,permDim));
scrambledim = trafo.scrambledim(1,permDim);
trafo.scrambledim = find(trafo.scrambledim(1,permDim) == 1);
transformDim = find(trafo.trafodim(2,permDim));
if(any(transformDim > testdim))
    transformDim = transformDim(transformDim <= testdim); % guarantess that e.g. a 3D recon is not applied to a 2D image
end
if(any(transformDim ~= 1:length(transformDim)))
    % permutation of dimensions is needed (place to be transformed
    % dimensions at the beginning), because not the first few dimensions
    % should be transformed
    tmp = 1:testdim;
    trafo.permRule = [transformDim, tmp(~ismember(tmp,transformDim))]; % permRule assumes dimensionalities of input kSpace: y-z-x-t-cha
    % catch case permRule = 1:testdim => no permutation needed
    if(all(trafo.permRule == tmp))
        trafo.permRule = [];
    else
        trafo.permRule = [trafo.permRule, max(trafo.permRule)+1]; % channels
    end
else
    trafo.permRule = [];
end
% trafo.transformDim = transformDim;
% N-dimensional transformation; 1D=1, 2D=2, 3D=3, 4D=4
trafo.shape = length(transformDim); 
if(length(transformDim) >= 5)
    msg = 'N-D, with N >= 5, sparsifying transformation is not supported';
    consExist = false;
    return
end
tmp = trafo.trafodim(1,permDim);
trafo.kspaceTrafoFFTdims = double(tmp);
trafo.kspaceTrafoFFTdims(trafo.kspaceTrafoFFTdims == 1) = find(tmp);
trafo.kspaceTrafoScrambledims = double(scrambledim);
trafo.kspaceTrafoScrambledims(trafo.kspaceTrafoScrambledims == 1) = find(scrambledim == 1);
if(~isempty(trafo.permRule))
    tmp = tmp(trafo.permRule(1:end-1));
    trafo.kspaceTrafoFFTdims = trafo.kspaceTrafoFFTdims(trafo.permRule(1:end-1));
    trafo.kspaceTrafoScrambledims = trafo.kspaceTrafoScrambledims(trafo.permRule(1:end-1));
end
trafo.kspaceTrafoFFTdims = trafo.kspaceTrafoFFTdims(1:trafo.shape);
trafo.kspaceTrafoFFTdims = trafo.kspaceTrafoFFTdims(trafo.kspaceTrafoFFTdims ~= 0);
trafo.kspaceTrafoScrambledims = trafo.kspaceTrafoScrambledims(1:trafo.shape);
trafo.kspaceTrafoScrambledims = trafo.kspaceTrafoScrambledims(trafo.kspaceTrafoScrambledims ~= 0);
if(isempty(trafo.fftdim) || (all(tmp(1:trafo.shape) == 0)) && ~isempty(trafo.kspaceTrafoFFTdims))
    trafo.kspaceTrafo = true;
else
    trafo.kspaceTrafo = false;
end
    
% set fftdim and scrambledim
if(~isempty(trafo.fftdim) && any(trafo.fftdim > testdim))
    trafo.fftdim = trafo.fftdim(trafo.fftdim <= testdim);
end
if(~isempty(trafo.scrambledim) && any(trafo.scrambledim > testdim))
    trafo.scrambledim = trafo.scrambledim(trafo.scrambledim <= testdim);
end
if(any(mod(kernelSize,2) == 0))
    msg = 'for kernelSize just odd sizes are allowed';
    consExist = false;
    return
end

% set reconDIM for Proximal averages
if(trafo.shape == 2 || trafo.shape == 1)
    reconDIM = '2D';
elseif(trafo.shape == 3)
    reconDIM = '3D';
end

% output cs_type, sparsifying transformation and along which dimensions
dimOutTmp = {'ky', 'kx', 'kz', 't'; 'y', 'x', 'z', 'f'};
dimStart = '';
dimEnd = '';
for i=1:testdim
    if(trafo.fftBA(1,i))
        dimStart = [dimStart, dimOutTmp{2,permDim(i)}, ' - '];
    else
        dimStart = [dimStart, dimOutTmp{1,permDim(i)}, ' - '];
    end
    if(trafo.trafodim(2,permDim(i)))
        dimEnd = [dimEnd, dimOutTmp{1+trafo.trafodim(1,permDim(i))+trafo.fftBA(1,i), permDim(i)}, 'B - '];
    else
        dimEnd = [dimEnd, dimOutTmp{1+trafo.trafodim(1,permDim(i))+trafo.fftBA(1,i), permDim(i)}, ' - '];
    end
end
fprintf('CS reconstruction: %s with %s, %s <=> %s\n', cstype, trafoType, dimStart(1:end-3), dimEnd(1:end-3));

if(~emptyPara(pwindow))
    fname = functions(window.type);
    fname = fname.function;        
    if(strcmp(fname,'rectwin'))
        window.windowingOn = false;
    else
        window.windowingOn = true;
    end
    switch fname
        case {'rectwin', 'barthannwin', 'bartlett', 'bohmanwin', 'parzenwin', 'triang'}
            % nop
        case {'blackman','blackmanharris', 'hamming', 'hann', 'flattopwin', 'nuttallwin', 'lanczos', 'welch', 'cosine'}
            if(~strcmp(window.windowOpt{1},'symmetric') && ~strcmp(window.windowOpt{1},'periodic'))
                msg = sprintf('Window option must be symmetric or periodic');
                consExist = false;
                return
            end
        case {'chebwin', 'tukeywin', 'gausswin', 'kaiser', 'planckTaper'}
            if(~strcmp(window.windowOpt{1},'symmetric') && ~strcmp(window.windowOpt{1},'periodic') || ~isscalar(window.windowOpt{2}) || window.windowOpt{2} <= 0)
                msg = sprintf('Window option must be symmetric or periodic and additional window parameter must be a scalar >0');
                consExist = false;
                return
            end
    end
else
    window.windowingOn = false;
end

if(~emptyPara(pespresso))
    if(allPara{1,pespresso}.state && allPara{1,pespresso}.iter < 0)
        msg = sprintf('ESPReSSo POCS iterations must be >0');
        consExist = false;
        return
    end
end

% consistency between recon type and chosen sparsifying transformation
% FOCUSS: all
% sparseMRI: FFT (p2DFT), Wavelet_lab
% SPIRiT, ESPIRiT: Wavelet_lab, FFT
% L1_Magic: FFT
% RecPF: FFT, Wavelet_lab
% Proximal: Wavelet_mat
if(strcmp(cstype,'FOCUSS'))
    if(~strcmp(trafoType,'fft') && ~strcmp(trafoType,'dct') && ~strcmp(trafoType,'wavelet_mat') && ~strcmp(trafoType,'pca') && ~strcmp(trafoType,'mellin') && ~strcmp(trafoType,'curvelet') && ~strcmp(trafoType,'surfacelet'))
        msg = sprintf('%s: not eligible transformation type for %s', trafoType, cstype);
        consExist = false;
        return
    end
elseif(strcmp(cstype,'sparseMRI'))
    if(~strcmp(trafoType,'fft') && ~strcmp(trafoType,'wavelet_lab') && ~strcmp(trafoType,'wavelet_mat'))
        msg = sprintf('%s: not eligible transformation type for %s', trafoType, cstype);
        consExist = false;
        return
    end
elseif(strcmp(cstype,'SPIRiT_CG') || strcmp(cstype,'SPIRiT_POCS') || strcmp(cstype,'ESPIRiT_CG') || strcmp(cstype,'ESPIRiT_L1'))
    if(~strcmp(trafoType,'fft') && ~strcmp(trafoType,'wavelet_lab'))
        msg = sprintf('%s: not eligible transformation type for %s', trafoType, cstype);
        consExist = false;
        return
    end
elseif(strcmp(cstype,'L1_Magic_TV') || strcmp(cstype,'L1_Magic_L1') || strcmp(cstype,'L1_Magic_L1Dantzig') || strcmp(cstype,'L1_Magic_TVDantzig'))
    if(~strcmp(trafoType,'fft'))
        msg = sprintf('%s: not eligible transformation type for %s', trafoType, cstype);
        consExist = false;
        return
    end
elseif(strcmp(cstype,'RecPF'))
    if(~strcmp(trafoType,'fft') && ~strcmp(trafoType,'wavelet_lab'))
        msg = sprintf('%s: not eligible transformation type for %s', trafoType, cstype);
        consExist = false;
        return
    end
elseif(strcmp(cstype,'FCSA_WaTMRI') || strcmp(cstype,'FCSA_SLEP') || strcmp(cstype,'FCSA_proxA') || strcmp(cstype,'BFCSA_proxA') || strcmp(cstype,'ADMM_proxA') || strcmp(cstype,'SB_proxA'))
    if(~strcmp(trafoType,'wavelet_mat'))
        msg = sprintf('%s: not eligible transformation type for %s', trafoType, cstype);
        consExist = false;
        return
    end   
end

if(~emptyPara(ptrafo))
    switch trafoType
        case 'fft'
            if(allPara{1,ptrafo}.windowing)
                fname = functions(allPara{1,ptrafo}.windowType);
                fname = fname.function;
                switch fname
                    case {'rectwin', 'barthannwin', 'bartlett', 'bohmanwin', 'parzenwin', 'triang'}
                        % nop
                    case {'blackman','blackmanharris', 'hamming', 'hann', 'flattopwin', 'nuttallwin', 'lanczos', 'welch', 'cosine'}
                        if(~strcmp(window.windowOpt{1},'symmetric') && ~strcmp(window.windowOpt{1},'periodic'))
                            msg = sprintf('Window option for fft must be symmetric or periodic');
                            consExist = false;
                            return
                        end
                    case {'chebwin', 'tukeywin', 'gausswin', 'kaiser', 'planckTaper'}
                        if(~strcmp(window.windowOpt{1},'symmetric') && ~strcmp(window.windowOpt{1},'periodic') || ~isscalar(window.windowOpt{2}) || window.windowOpt{2} <= 0)
                            msg = sprintf('Window option for fft must be symmetric or periodic and additional window parameter must be a scalar >0');
                            consExist = false;
                            return
                        end
                    otherwise
                        msg = sprintf('Unknown window type for fft');
                        consExist = false;
                        return
                end
                if(trafo.shape == 3)
                    msg = '3D windowing not yet supported';
                    consExist = false;
                    return
                end
            end  
            
        case 'dct'
            if(strcmp(sPrecision,'single'))
                msgWarn = [msgWarn,'Switching calculation precision: double precision needed for dct; '];
                sPrecision = 'double';
            end
            
        case 'wavelet_mat'
            if(~strcmp(allPara{1,ptrafo}.waveletFilter,'haar') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'db') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'coif') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'sym') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'dmey') && ...
                    ~strcmp(allPara{1,ptrafo}.waveletFilter,'bior') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'rbio'))
                msg = 'No valid Wavelet QMF filter chosen';
                consExist = false;
                return
            end
            switch allPara{1,ptrafo}.waveletFilter
                case 'haar'
                    if(allPara{1,ptrafo}.waveletFilterSize ~= 1)
                        msg = 'Wavelet filter size must be 1';
                        consExist = false;
                        return
                    end
                case 'db'
                    if(~any(ismember(1:45,allPara{1,ptrafo}.waveletFilterSize)))
                        msg = 'Wavelet filter size must be between 1 and 45';
                        consExist = false;
                        return
                    end
                case 'coif'
                    if(~any(ismember(1:5,allPara{1,ptrafo}.waveletFilterSize)))
                        msg = 'Wavelet filter size must be between 1 and 5';
                        consExist = false;
                        return
                    end
                case 'sym'
                    if(~any(ismember(2:45,allPara{1,ptrafo}.waveletFilterSize)))
                        msg = 'Wavelet filter size must be between 2 and 45';
                        consExist = false;
                        return
                    end
                case 'dmey'
                    if(allPara{1,ptrafo}.waveletFilterSize <= 0)
                        msg = 'Wavelet filter size must be larger than 0';
                        consExist = false;
                        return
                    end
                case {'bior', 'rbior'}
                    if(~any(ismember([1.1, 1.3, 1.5, 2.2, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.7, 3.9, 4.4, 5.5, 6.8],allPara{1,ptrafo}.waveletFilterSize)))
                        msg = 'Wavelet filter size not supported';
                        consExist = false;
                        return
                    end
            end
            if(~strcmp(allPara{1,ptrafo}.extMode,'zpd') && ~strcmp(allPara{1,ptrafo}.extMode,'sym') && ~strcmp(allPara{1,ptrafo}.extMode,'symh') && ~strcmp(allPara{1,ptrafo}.extMode,'symw') && ...
                    ~strcmp(allPara{1,ptrafo}.extMode,'asym') && ~strcmp(allPara{1,ptrafo}.extMode,'asymh') && ~strcmp(allPara{1,ptrafo}.extMode,'asymw') && ~strcmp(allPara{1,ptrafo}.extMode,'spd') && ...
                    ~strcmp(allPara{1,ptrafo}.extMode,'sp1') && ~strcmp(allPara{1,ptrafo}.extMode,'sp0') && ~strcmp(allPara{1,ptrafo}.extMode,'ppd'))
                msg = 'Unknown extension mode';
                consExist = false;
                return
            end
            if(allPara{1,ptrafo}.waveletStages <= 0)
                msg = 'Amount of decomposition stages must be larger than 0';
                consExist = false;
                return
            end
            if(trafo.shape == 4)
                msg = '4D sparsifying transformation not supported yet!';
                consExist = false;
                return
            end
            
        case 'wavelet_lab'
            if(trafo.shape == 1 || trafo.shape == 3 || trafo.shape == 4)
               msg = sprintf('Full %dD transformation not supported by Wavelab850', trafo.shape);
               consExist = false;
               return
            end
            if(strcmp(dimension,'3D') || strcmp(dimension,'4D') || (strcmp(dimension,'2D') && dim(4) > 1)) 
                msgWarn = [msgWarn,'Wavelab850 just supports 2D images. Iteration over remaining dimensions.; '];
            end
            if(~strcmp(allPara{1,ptrafo}.waveletFilter,'Haar') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'Beylkin') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'Coiflet') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'Daubechies') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'Symmlet') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'Vaidyanathan') && ~strcmp(allPara{1,ptrafo}.waveletFilter,'Battle'))
                msg = 'No valid Wavelet QMF filter chosen';
                consExist = false;
                return
            end
            
            if((strcmp(cstype,'SPIRiT_CG') || strcmp(cstype,'SPIRiT_POCS') || strcmp(cstype,'ESPIRiT_CG') || strcmp(cstype,'ESPIRiT_L1')) && (~allPara{1,ptrafo}.flagThresholding))
                trafo.flagThresholding = true;
                msgWarn = [msgWarn,'For SPIRiT/ESPIRiT wavelet thresholding must be enabled!; '];
            end
                      
            [allPara{1,ptrafo}.waveletFilter,TI]=strtok(allPara{1,ptrafo}.waveletFilter,'_'); % cut off TI-flag
            
            if(~isempty(TI) && ~strcmp(TI,'TI'))
                msg = 'Wavelet lab just allows ''TI'' flag';
                consExist = false;
                return
            end
            
            switch allPara{1,ptrafo}.waveletFilter
                case {'Haar', 'Beylkin', 'Vaidyanathan'}
                    if(allPara{1,ptrafo}.waveletFilterSize <=0)
                        msg = 'Wavelet filter size must be larger than 0';
                        consExist = false;
                        return
                    end
                case 'Coiflet'
                    if(~any(ismember(1:5,allPara{1,ptrafo}.waveletFilterSize)))
                        msg = 'Wavelet filter size must be between 1 and 5';
                        consExist = false;
                        return
                    end
                case 'Daubechies'
                    if(~any(ismember(4:2:20,allPara{1,ptrafo}.waveletFilterSize)))
                        msg = 'Wavelet filter size must be 4, 6, 8, 10, 12, 14, 16, 18 or 20';
                        consExist = false;
                        return
                    end
                case 'Symmlets'
                    if(~any(ismember(4:10,allPara{1,ptrafo}.waveletFilterSize)))
                        msg = 'Wavelet filter size must be between 4 and 10';
                        consExist = false;
                        return
                    end
                case 'Battle'
                    if(~any(ismember([1,3,5],allPara{1,ptrafo}.waveletFilterSize)))
                        msg = 'Wavelet filter size must be 1, 3 or 5';
                        consExist = false;
                        return
                    end
            end
    
            if(allPara{1,ptrafo}.waveletStages <= 0)
                msg = 'Wavelet decomposition level should be greater than 0';
                consExist = false;
                return
            end
            
            if(allPara{1,ptrafo}.wavWeight <= 0)
                msg = 'Wavelet soft threshold should be greater than 0';
                consExist = false;
                return
            end
            
        case 'mellin'
            if(trafo.shape ~= 2 && trafo.shape ~= 3)
                msg = sprintf('%s reconstruction not supported', trafo.shape);
                consExist = false;
                return
            end
            if(isinf(allPara{1,ptrafo}.extrapVal) || isnan(allPara{1,ptrafo}.extrapVal))
                msg = 'Interpolation value must be a real-valued number';
                consExist = false;
                return
            end
            if(any(allPara{1,ptrafo}.trafoSize <= 0) || ~all(mod(allPara{1,ptrafo}.trafoSize,1) == 0))
                msg = 'Trafo size must be larger than 0 and an integer multiple of the original dimension';
                consExist = false;
                return
            end
            if(~strcmp(allPara{1,ptrafo}.interp,'nearest') && ~strcmp(allPara{1,ptrafo}.interp,'linear') && ~strcmp(allPara{1,ptrafo}.interp,'spline'))
                msg = 'Interpolation type must be nearest or linear';
                consExist = false;
                return
            end
            if(allPara{1,ptrafo}.sigma <= 0)
                msg = 'Distance scaling must be larger than 0';
                consExist = false;
                return
            end
            
        case 'curvelet' 
            if(trafo.shape ~= 2 && trafo.shape ~= 3)
                msg = sprintf('%s reconstruction not supported', trafo.shape);
                consExist = false;
                return
            end
            if(~(strcmp(dimension,'2D') && all(trafo.trafodim(2,1:2))) && ...
                    ~((strcmp(dimension,'3D') || strcmp(dimension,'4D')) && trafo.shape == 2 && all(trafo.trafodim(2,1:2))) && ...
                    ~((strcmp(dimension,'3D') || strcmp(dimension,'4D')) && trafo.shape == 3 && all(trafo.trafodim(2,1:3))))
                msgWarn = [msgWarn,'curvelet is better suited for aligned image dimensions; '];
            end
            if(allPara{1,poversampling}{2,1} == 0)
                msgWarn = [msgWarn,'Anti-Aliasing correction along x dimension should be performed before the reconstruction in order to avoid large zero-padding; '];
            end
            if (allPara{1,ptrafo}.allCurvelets ~=0 && allPara{1,ptrafo}.finest ~=1)
                msg = 'realVal has to be either 0 curvelets or 1 for wavelets at minimal scale';
                consExist = false;
                return
            end
            if (mod(allPara{1,ptrafo}.nbdstz_coarse, 4) > 0)
                msg = 'angular spacing has to be an integer multiple of 4';
                consExist = false;
                return
            end
            if(ispc)
                msg = 'curvelet is not supported for Windows machines';
                consExist = false;
                return
            end           
            
        case 'surfacelet'
            if(trafo.shape ~= 2 && trafo.shape ~= 3)
                msg = sprintf('%dD reconstruction not supported', trafo.shape);
                consExist = false;
                return
            end
            if(~(strcmp(dimension,'2D') && all(trafo.trafodim(2,1:2))) && ...
                    ~((strcmp(dimension,'3D') || strcmp(dimension,'4D')) && trafo.shape == 2 && all(trafo.trafodim(2,1:2))) && ...
                    ~((strcmp(dimension,'3D') || strcmp(dimension,'4D')) && trafo.shape == 3 && all(trafo.trafodim(2,1:3))))
                msgWarn = [msgWarn,'surfacelet is better suited for aligned image dimensions; '];
            end
            if(~allPara{1,poversampling}{2,1})
                msgWarn = [msgWarn,'Anti-Aliasing correction along x dimension should be performed before the reconstruction in order to avoid large zero-padding; '];
            end
            if(~any(ismember([1, 1.5, 2],allPara{1,ptrafo}.Pyr_mode)))
                msg = 'Pyr_mode must be 1, 1.5 or 2';
                consExist = false;
                return
            end
            if(~strcmp(allPara{1,ptrafo}.HGfname,'ritf'))
                msg = 'Unknown hourglass filter bank';
                consExist = false;
                return
            end
            if(allPara{1,ptrafo}.bo <= 0 || allPara{1,ptrafo}.msize <= 0 || allPara{1,ptrafo}.beta <= 0 || allPara{1,ptrafo}.lambda <= 0 || allPara{1,ptrafo}.decompLevels <= 0)
                msg = 'Transformation parameter must be larger than 0';
                consExist = false;
                return
            end   
            % correct dimensionality of Lev_array
            if((size(allPara{1,ptrafo}.Lev_array{1},1) < 2 && trafo.shape == 2) || (size(allPara{1,ptrafo}.Lev_array{1},1) < 3 && trafo.shape == 3))
                msg = 'Lev_array dimensionality is set too small';
                consExist = false;
                return
            end
            if(size(allPara{1,ptrafo}.Lev_array{1},1) > 2 && trafo.shape == 2) 
                trafo.Lev_array = cellfun(@(x) x(1:2,1:2), trafo.Lev_array, 'UniformOutput', false);
                msgWarn = [msgWarn,'Lev_array dimensionality is set too large. Reducing to 2D; '];
            elseif((size(allPara{1,ptrafo}.Lev_array{1},1) > 3 && trafo.shape == 3))
                trafo.Lev_array = cellfun(@(x) x(1:3,1:3), trafo.Lev_array, 'UniformOutput', false);
                msgWarn = [msgWarn,'Lev_array dimensionality is set too large. Reducing to 3D; '];
            end
            % correct msize if decomposed image is getting too small
            downSamp = cell(1,trafo.decompLevels+1); % cell 1: inital downsampling due to Lev_array{1}, cell 2: downsampling from stage 1 to stage 2 due to pyr_mode and Lev_array{2}
            [downSamp{:}] = deal(zeros(size(trafo.Lev_array{1})));
            [downSamp{3:trafo.decompLevels+1}] = deal(ones(size(trafo.Lev_array{1})));
            for i=1:trafo.decompLevels
                lidx = trafo.Lev_array{i} ~= 0 & trafo.Lev_array{i} ~= -1;
                if(any(lidx))
                    downSamp{i}(lidx) = downSamp{i}(lidx) + trafo.Lev_array{i}(lidx);
                end
            end
            downSamp = cell2mat(shiftdim(cellfun(@(x) 2.^x, downSamp, 'UniformOutput', false),-1));
            downSamp(:,:,2) = 2.^(log2(downSamp(:,:,2)) + log2(trafo.Pyr_mode));
            trafo.downSamp = max(prod(downSamp,3),[],1); 
            trafo.downSamp = trafo.downSamp(1:trafo.shape);
%             if(trafo.Pyr_mode == 1.5)
%                 trafo.downSamp = lcm(trafo.downSamp,3); % from source code of surfaclet decomposition (for Lev_array{0} => image size must be dividable by 2^Lev_array and 3)
%             end
%             if(trafo.Pyr_mode == 2)
%                 trafo.downSamp = 2*trafo.downSamp;
%             end
            % find minimal padsize
            maptmp = [2 1 3 4]; % mapping x-y-z(-t) to y-x-z(-t)
            if(~isempty(trafo.permRule))
                minsize = permute(dim,trafo.permRule(1:end-1));
                maptmp  = permute(maptmp,trafo.permRule(1:end-1));
            else
                minsize = dim; % y-x-z-t-cha
            end
            for i=1:length(maptmp)
                if(maptmp(i) == 4)
                    continue;
                end
                if(allPara{1,poversampling}{2,i})
                    minsize(maptmp(i)) = numel(allPara{1,poversampling}{1,i}); % x-y-z
                end
            end
            minsize = minsize(1:trafo.shape);
%             if(~isempty(allPara{1,poversampling}{1,2})), minsize(1) = minsize(1) / numel(allPara{1,poversampling}{1,2}); end;
%             if(~isempty(allPara{1,poversampling}{1,1})), minsize(2) = minsize(2) / numel(allPara{1,poversampling}{1,1}); end;
%             if(~isempty(allPara{1,poversampling}{1,3})), minsize(3) = minsize(3) / numel(allPara{1,poversampling}{1,3}); end;
            minsize(minsize == 1) = [];
            n = minsize./trafo.downSamp;
            if(any(n < 1))
                n(n<1) = 1;
                msgWarn = [msgWarn,'Surfacelet parameters are set too strict -> large zero padding. Consider rescaling!; '];
            end
            n = ceil(n);
            trafo.minsize = max(n.*trafo.downSamp);
            minsize = trafo.minsize/max(trafo.downSamp);
            if(minsize < 2*trafo.msize - 1)
                trafo.msize = floor((minsize + 1)/2);
            end           
            if(~allPara{1,pflagZeropadding})
                trafo.padOnSameSize = false;
            end
    end
end
    
if((strcmp(cstype,'SPIRiT_CG') || strcmp(cstype,'SPIRiT_POCS') || strcmp(cstype,'ESPIRiT_CG') || strcmp(cstype,'ESPIRiT_L1')) && ~strcmp(trafoType,'wavelet_lab'))
    msg = 'SPIRiT/ESPIRiT algorithm just works with set parameters for Wavelab850';
    consExist = false;
    return
end

if((strcmp(cstype,'SPIRiT_CG') || strcmp(cstype,'SPIRiT_POCS') || strcmp(cstype,'ESPIRiT_CG') || strcmp(cstype,'ESPIRiT_L1')) && strcmp(sPrecision,'single'))
    msgWarn = [msgWarn,sprintf('Switching calculation precision: double precision needed for %s; ',cstype)];
    sPrecision = 'double';
end

if(strcmp(cstype,'FCSA_TV') || strcmp(cstype,'FCSA_NLTV'))
    msg = 'Deprecated code: No longer supported';
    consExist = false;
    return
end

if((strcmp(cstype,'FCSA_TV') || strcmp(cstype,'FCSA_NLTV')) && ~strcmp(trafoType,'wavelet_mat'))
    msg = 'FCSA algorithm just works with set parameters for MatLab Wavelet toolbox';
    consExist = false;
    return
end

% if((strcmp(cstype,'ESPIRiT_CG') || strcmp(cstype,'ESPIRiT_L1')) && ~ispc && ~allPara{1,pflagToolbox})
%     msgWarn = [msgWarn,'ESPIRiT on linux machine should use the ESPIRiT Linux toolbox; '];
% elseif((strcmp(cstype,'ESPIRiT_CG') || strcmp(cstype,'ESPIRiT_L1')) && ispc && allPara{1,pflagToolbox})
%     flagToolbox = false;
%     msgWarn = [msgWarn,'ESPIRiT on windows machine cannot use the ESPIRiT Linux toolbox; '];
% end

if((strcmp(cstype,'L1_Magic_TV') || strcmp(cstype,'L1_Magic_L1') || strcmp(cstype,'L1_Magic_L1Dantzig') || ...
        strcmp(cstype,'L1_Magic_TVDantzig')) && ~strcmp(trafoType,'fft'))
    msg = sprintf('%s currently just works with fourier transformation',cstype);
    consExist = false;
    return
end

% adapt flags
if(allPara{1,pLambdaGroup} > 0)
    flags.flagGroup = true;
else
    flags.flagGroup = false;
end
if(allPara{1,pLambdaNLTV} > 0)   
    flags.flagNTLV = true;
else
    flags.flagNTLV = false;
end
if(allPara{1,pLambdaTV} > 0)
    flags.flagTV = flags.flagTV & true;
    flags.flagNLTV = flags.flagNLTV & true;
else
    flags.flagTV = false;
end

if(strcmp(cstype,'FCSA_WaTMRI') || strcmp(cstype,'FCSA_SLEP') || strcmp(cstype,'FCSA_proxA') || strcmp(cstype,'BFCSA_proxA') || strcmp(cstype,'ADMM_proxA') || strcmp(cstype,'SB_proxA'))
    if(trafo.shape == 2) 
        reconDIM = '2D';
    elseif(trafo.shape == 3)
        reconDIM = '3D';
        if(flags.flagTV_iso)
            flags.flagTV_iso = false;
            msgWarn = [msgWarn, 'Switching isotropic TV to anisotropic; '];
        end
        if(flags.flagGroup)
            msgWarn = [msgWarn, 'Group thresholding not fully supported for 3D. Looping over third dimension; '];
        end
    else
        msg = sprintf('Reconstruction dimensionality %dD for algorithm %s not supported', trafo.shape, cstype);
        consExist = false;
        return
    end
end

if(strcmp(cstype,'FCSA_WaTMRI'))
    if(strcmp(reconDIM,'3D'))
        msgWarn = [msgWarn, 'Switching reconstruction dimensionality to 2D; '];
        reconDIM = '2D';
    end
    if(~flags.flagRealOnly)
        msgWarn = [msgWarn, 'Switching to real-valued calculations; '];
        flags.flagRealOnly = true;
    end
elseif(strcmp(cstype,'FCSA_SLEP'))
    if(strcmp(reconDIM,'3D'))
        msgWarn = [msgWarn, 'Switching reconstruction dimensionality to 2D; '];
        reconDIM = '2D';
    end
    if(~flags.flagRealOnly)
        msgWarn = [msgWarn, 'Switching to real-valued calculations; '];
        flags.flagRealOnly = true;
    end
end               


end
