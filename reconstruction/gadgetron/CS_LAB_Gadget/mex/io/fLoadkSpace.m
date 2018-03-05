function [dKSpace, measPara, file] = fLoadkSpace(path)
%% Load kSpace data
%
% (c) Martin Schwartz and Thomas Kuestner
% -------------------------------------------------------------------------

    file='';
    warning('off','MATLAB:dispatcher:nameConflict');
    para = 0;
    % prepare m-files[
    currpath = fileparts(mfilename('fullpath'));

    % load ADCs
    if(ischar(path))
        if(isdir(path))
            [file, path] = uigetfile({'*.mat;*.dat', 'ADC measdata (*.mat, *.dat)'} ,'Select ADC measdata file', 'MultiSelect', 'off', path);
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
   
    % read inputs
    if(lReadFromFile)
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
    
    %% get data dimension and loop counters
    espresso.reconType = 1;
    postproc.aniso = [];
    [measPara, ~, ~] = extractMeasPara( iLC, measPara, espresso, postproc, logical([0 0 0]), drecksMDH, iLCPositions );

    %% extract kSpace data from measurement file
    [kSpace, ~] = extractKSpace_main(dData, iLC, iEvalInfoMask, measPara);
    
    dKSpace = zeros([size(kSpace{1}),size(kSpace,2)]);
    for iI = 1:size(kSpace,2)
        dKSpace(:,:,:,iI) = kSpace{iI};
    end
    % y-x-z-c -> x-y-z-c
    dKSpace = permute(dKSpace, [2 1 3 4]);
    
end