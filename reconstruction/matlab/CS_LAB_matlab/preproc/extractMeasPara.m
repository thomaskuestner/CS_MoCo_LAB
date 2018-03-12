function [ measPara, espresso, postproc ] = extractMeasPara( iLC, measPara, espresso, postproc, flagOversampling, drecksMDH, iLCPositions )
% extract measurement parameter from the loop counter and from the
% DrecksMDH (20 samples ADC)

% output:
% measPara      struct containing all measurement parameters (including
%               dimensions)
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------


%% dim and LCall
% dim: nPha, nFreq, nZ, nTime, nCha
% LCall: nSlices, nAcquisitions, nEchos, nRepetitions
if(~isempty(iLC))
    measPara.dim = ones(1,5);
    measPara.LCall = [0 1 1 1];

    measPara.dim(1) = double(max(iLC(:,3)) - min(iLC(:,3)) + 1); % nPha    = 'Lin'
    measPara.dim(2) = double(iLC(1,1));                          % nFreq   = 'Smp'
    measPara.dim(3) = double(max(iLC(:,6)) - min(iLC(:,6)) + 1); % nZ      = 'Par'
    measPara.dim(4) = double(max(iLC(:,8)) - min(iLC(:,8)) + 1); % nTime   = 'Pha' -- ATTENTION: newer FLASH implentation: measurements = time => CORRECT IT!!!
    measPara.dim(5) = double(iLC(1,2));                          % nCha    = 'Cha'

    measPara.LCall(1) = double(max(iLC(:,5)) - min(iLC(:,5)) + 1); % nSlices       = 'Sli'
    measPara.LCall(2) = double(max(iLC(:,4)) - min(iLC(:,4)) + 1); % nAcquisitions = 'Acq' / nAverages = 'Avg'
    measPara.LCall(3) = double(max(iLC(:,7)) - min(iLC(:,7)) + 1); % nEchos        = 'Eco' % TODO: check if working correctly
    measPara.LCall(4) = double(max(iLC(:,9)) - min(iLC(:,9)) + 1); % nRepetitions  = 'Rep'

    if(measPara.dim(3) == 1)
        measPara.dimension = '2D';
    else
        if(measPara.dim(4) == 1)
            measPara.dimension = '3D';
        else
            measPara.dimension = '4D';
        end    
    end
end

%% ESPReSSo
if(~isempty(iLC) && max(iLC(:,14)) > 81) % see IDEA: WIP_UI.h
    espresso.state = true;
    pfn_mapping = [90, 99, 108, 117; 0.5, 0.625, 0.75, 0.875];
    espresso.pfn = pfn_mapping(2,pfn_mapping(1,:) == max(iLC(:,14)));
    % correct dimensionality
    measPara.dim(1) = double(max(iLC(:,3)) + 1);
    measPara.dim(3) = double(max(iLC(:,6)) + 1);
else
    if(~isfield(espresso,'state'))
        espresso.state = false;
    end
    if(~isfield(espresso,'pfn'))
        espresso.pfn = 0;
    end
    if(~isfield(espresso,'direction'))
        espresso.direction = 'off';
    end
end


%% oversampling
% initialize correction for oversampling along x-y-z/slice direction
if(~isfield(measPara,'oversampling'))
    measPara.oversampling = cell(2,3);
    measPara.oversampling(2,:) = num2cell(flagOversampling);
    if(isfield(measPara,'dim'))
        % frequency oversampling
        measPara.oversampling{1,1} = measPara.dim(2)/4:measPara.dim(2)*3/4 - 1; % set cut-out x-region for anti-aliasing
        % phase oversampling
        if(~isempty(iLC))
            phaseOver = double(max(iLC(:,13)))/1000;
        else
            phaseOver = 0;
        end
        % phaseDiff = abs(measPara.dim(1) - round(fzero(@(x) x+round(phaseOver*x)-measPara.dim(1),measPara.dim(1))));
        phaseDiff = abs(measPara.dim(1) - round(fzero(@(x) round(x*(1+phaseOver)) - measPara.dim(1), measPara.dim(1))));
        if(measPara.oversampling{2,2})
            measPara.oversampling{1,2} = round(phaseDiff/2)+1:(measPara.dim(1)-round(phaseDiff/2));
        else
            measPara.oversampling{1,2} = 1:measPara.dim(1);
        end
        % slice oversampling
        if(~isempty(iLC))
            sliceOver = double(max(iLC(:,12)))/1000;
        else
            sliceOver = 0;
        end
        % if(strcmp(measPara.dimension,'2D'))
        %     sliceDiff = abs(measPara.LCall(1) - round(fzero(@(x) x+round(sliceOver*x)-measPara.LCall(1),measPara.LCall(1))));
        %     measPara.oversampling{1,3} = round(sliceDiff/2)+1:(measPara.LCall(1)-round(sliceDiff/2));
        % else
        if(strcmp(measPara.dimension,'3D') || strcmp(measPara.dimension,'4D')) % no slice oversampling for 2D
        %     sliceDiff = abs(measPara.dim(3) - round(fzero(@(x) x+round(sliceOver*x)-measPara.dim(3),measPara.dim(3))));
            sliceDiff = abs(measPara.dim(3) - round(fzero(@(x) round(x*(1+sliceOver))-measPara.dim(3),measPara.dim(3))));
            if(measPara.oversampling{2,3})
                measPara.oversampling{1,3} = round(sliceDiff/2)+1:(measPara.dim(3)-round(sliceDiff/2));
            else
                measPara.oversampling{1,3} = 1:measPara.dim(3);
            end
        else
            measPara.oversampling{1,3} = 1:measPara.LCall(1);
        end
    end
end

%% anisotropic correction
% not possible without additional information
measPara.aniso = [];


%% extract DrecksMDH parameters
if(~isempty(drecksMDH))
    if(~isstruct(drecksMDH)) % drecksMDH 1.0
        measPara = fMeasExtractMeasPara( drecksMDH, measPara, iLC, iLCPositions );
    else % drecksMDH 2.0 -> new drecksMDH style -> convert to measPara-struct
        measPara = convertToMeasPara(drecksMDH, measPara); % measPara now just contains all really necessary variables
    end
    espresso.direction = measPara.ESPReSSoDirection;
    espresso.pfn    = measPara.ESPReSSo;    
end % else no drecksMDH at all

% ESPReSSo
if(isempty(espresso.pfn) || espresso.pfn == 1 || strcmp(espresso.direction,'off'))
    espresso.state = false;
else
    espresso.state = true;
end

if(~isempty(postproc.aniso))
    warning('Anisotropic correction is set to manual');
    tmpAniso = postproc.aniso([2 1 3]);
    tmpUnequal = cellfun(@(x) length(x), measPara.oversampling(1,:)) ~= tmpAniso;
    if(any(tmpUnequal))
        for i=1:find(tmpUnequal)
            if(i==1) % no oversampling correction for x direction
                continue;
            end
            measPara.oversampling{1,i} = 1:tmpAniso;
        end
    end
else
    postproc.aniso = measPara.aniso;
end


end
