function measPara = convertToMeasPara( drecksMDH, measPara )
%CONVERTTOMEASPARA convert from drecksMDH 2.0-struct to measPara-struct,
% but measPara-struct now contains just all needed variables
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

if(drecksMDH.Seq.Is3D)
    measPara.dimensionSeq = '3D';
else
    measPara.dimensionSeq = '2D';
end
measPara.sliceThickness = drecksMDH.Geo.FOV(3)/drecksMDH.Geo.MatrixSize(3);
measPara.lines = drecksMDH.Geo.MatrixSize(2);
measPara.partitions = drecksMDH.Geo.MatrixSize(3);
measPara.baseResolution = drecksMDH.Geo.MatrixSize(1);
measPara.FOVread = drecksMDH.Geo.FOV(1);
measPara.FOVphase = drecksMDH.Geo.FOV(2);
measPara.sequenceName = drecksMDH.Seq.Sequence;
if(isfield(drecksMDH.Geo,'FFTLength'))
    measPara.ftlenPhase = drecksMDH.Geo.FFTLength(1);
    measPara.ftlenPar = drecksMDH.Geo.FFTLength(2);
else
    measPara.ftlenPhase = [];
    measPara.ftlenPar = [];
end
measPara.imageWidth = [];
measPara.imageHeight = [];

if(isfield(drecksMDH,'Accel') && isfield(drecksMDH.Accel,'EspressoDir'))
    dir_mapping = {'y', 'z', 'off'};
    measPara.ESPReSSoDirection = dir_mapping{drecksMDH.Accel.EspressoDir};
    pfn_mapping = [81, 90, 99, 108, 117; 1, 0.5, 0.625, 0.75, 0.875];
    measPara.ESPReSSo = pfn_mapping(2,pfn_mapping(1,:) == drecksMDH.Accel.EspressoFactor);    
else
    measPara.ESPReSSoDirection = 'off';
    measPara.ESPReSSo = 1;
end

if(isfield(drecksMDH,'Wip'))
    if(isfield(drecksMDH.Wip,'SamplingFactor'))
        measPara.CSAcceleration = drecksMDH.Wip.SamplingFactor;
    end
    if(isfield(drecksMDH.Wip,'Phases'))
       	measPara.dim(4) = drecksMDH.Wip.Phases;
    end    
    if(isfield(drecksMDH.Wip,'FullySampled'))
        measPara.CSFullySampled = drecksMDH.Wip.FullySampled;
    end
    if(isfield(drecksMDH.Wip,'Mask'))
        switch drecksMDH.Wip.Mask(1:2)
            case '06', measPara.CSSamplingType = 'Poisson';
            case '07', measPara.CSSamplingType = 'Random';
            case '08', measPara.CSSamplingType = 'Gaussian';
        end
        if(length(drecksMDH.Wip.Mask) > 3)
            switch drecksMDH.Wip.Mask(4:5)
                case '01', measPara.CSVDmap = 'None';
                case '02', measPara.CSVDmap = 'Point';
                case '03', measPara.CSVDmap = 'Block';
                case '04', measPara.CSVDmap = 'Ellipse';
                case '05', measPara.CSVDmap = 'Ring';
            end
    %         switch drecksMDH.Wip.Mask(7:8)
    %             case '14', measPara.CSBodyRegion = 'None';
    %             case '15', measPara.CSBodyRegion = 'Head';
    %             case '16', measPara.CSBodyRegion = 'Thorax';
    %             case '17', measPara.CSBodyRegion = 'Abdomen';
    %             case '18', measPara.CSBodyRegion = 'Pelvis';
    %         end
        end
    end
end

if(isfield(drecksMDH.Seq,'Bandwidth'))
    measPara.bandwidthPerPixel = drecksMDH.Seq.Bandwidth;
end

if(isfield(drecksMDH.Geo,'NormSag'))
    if(drecksMDH.Geo.NormSag)
        measPara.orientation    = 'sagittal';
    elseif(drecksMDH.Geo.NormCor)
        measPara.orientation    = 'coronal';
    elseif(drecksMDH.Geo.NormTra)
        measPara.orientation    = 'transversal';
    else
        measPara.orientation = 'rotated';
    end

    % pixel positions
    if(isfield(drecksMDH.Geo,'Shift'))
        measPara.position = drecksMDH.Geo.Shift;
    else
        measPara.position = zeros(1,3);
    end
    measPara.pixelSpacing = [drecksMDH.Geo.FOV(2)/drecksMDH.Geo.MatrixSize(2), drecksMDH.Geo.FOV(1)/drecksMDH.Geo.MatrixSize(1), drecksMDH.Geo.FOV(3)/drecksMDH.Geo.MatrixSize(3)]; % y-x-z/slice
    center = [drecksMDH.Geo.MatrixSize(2)/2, drecksMDH.Geo.MatrixSize(1)/2, drecksMDH.Geo.MatrixSize(3)/2];

    switch measPara.orientation
        case 'sagittal'
            measPara.xPos = measPara.position(1) + ([1:drecksMDH.Geo.MatrixSize(3)].' - center(3)) .* measPara.pixelSpacing(3);
            measPara.yPos = measPara.position(2) + ([1:drecksMDH.Geo.MatrixSize(2)].' - center(1)) .* measPara.pixelSpacing(1);
            measPara.zPos = measPara.position(3) + ([1:drecksMDH.Geo.MatrixSize(1)].' - center(2)) .* measPara.pixelSpacing(2);
        case 'coronal'
            measPara.xPos = measPara.position(1) + ([1:drecksMDH.Geo.MatrixSize(1)].' - center(2)) .* measPara.pixelSpacing(2);
            measPara.yPos = measPara.position(2) + ([1:drecksMDH.Geo.MatrixSize(3)].' - center(3)) .* measPara.pixelSpacing(3);
            measPara.zPos = measPara.position(3) + ([1:drecksMDH.Geo.MatrixSize(2)].' - center(1)) .* measPara.pixelSpacing(1);
        case 'transversal'
            measPara.xPos = measPara.position(1) + ([1:drecksMDH.Geo.MatrixSize(1)].' - center(2)) .* measPara.pixelSpacing(2);
            measPara.yPos = measPara.position(2) + ([1:drecksMDH.Geo.MatrixSize(2)].' - center(1)) .* measPara.pixelSpacing(1);
            measPara.zPos = measPara.position(3) + ([1:drecksMDH.Geo.MatrixSize(3)].' - center(3)) .* measPara.pixelSpacing(3);
        otherwise 
            %
    end
end

% correct values
if(~isfield(measPara,'dim'))
    measPara.dim = zeros(1,5);
    measPara.LCall = zeros(1,4);
end
measPara.dim(1) = measPara.lines;
measPara.dim(2) = drecksMDH.Geo.OverSampling(1) * measPara.baseResolution; % due to CS_Trufi prescans
measPara.dim(3) = measPara.partitions;
% if(strcmp(measPara.sequenceName, 'CS_Trufi') && ~isempty(iLC)) % number of kernels (SetLC) = cardiac phases = time domain
%     measPara.dim(4) = double(max(iLC(:,10)) - min(iLC(:,10)) + 1); % SetLC
% end
measPara.LCall(2) = drecksMDH.LC.Averages;
if(measPara.dim(3) <= 1)
    measPara.dimension = '2D';
else
    if(measPara.dim(4) > 1)
        measPara.dimension = '4D';
    else
        measPara.dimension = '3D';
    end    
end
measPara.dim(6) = 1;
if(isfield(drecksMDH,'Wip'))
    if(isfield(drecksMDH.Wip,'NSamples') && drecksMDH.Wip.NSamples > 0)
        measPara.dimension = '5D';
        measPara.dim(6) = drecksMDH.Wip.NSamples;
    elseif(isfield(drecksMDH.Wip,'NCardGates') && drecksMDH.Wip.NCardGates > 0)
        measPara.dimension = '5D';
        measPara.dim(6) = drecksMDH.Wip.NCardGates;        
    end
end
if(isfield(measPara, 'sequenceName') && strcmp(measPara.sequenceName,'CS_Retro'))
    measPara.CSAcceleration = 2; % something unequal to 1
end

% local variables, just needed for further calculations
FOVphasePercentage = measPara.FOVphase/measPara.FOVread;
phaseOversampling = drecksMDH.Geo.OverSampling(2) - 1;
sliceOversampling = drecksMDH.Geo.OverSampling(3) - 1;

%% oversampling and anisotropy

% OVERSAMPLING PHASE INCORRECT!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% initialize anisotropic correction (y-x-z/slice-t)
measPara.aniso = [measPara.dim(1), measPara.dim(2), measPara.dim(3), measPara.dim(4), measPara.dim(6)];
% oversampling correction exactly possible => no case switch necessary later
if(~isfield(measPara,'oversampling'))
    measPara.oversampling = cell(2,3);
end
measPara.oversampling{2,2} = true;
measPara.oversampling{2,3} = true;

% readout oversampling
readoutDiff = (drecksMDH.Geo.OverSampling(1) - 1) * measPara.baseResolution;
measPara.oversampling{1,1} = readoutDiff/2+1:measPara.dim(2)-readoutDiff/2;
if(isempty(measPara.imageWidth) || measPara.imageWidth == 0)
    measPara.imageWidth = length(measPara.oversampling{1,1});
end

% correct anisotropy if oversampling correction was done during recon
if(measPara.oversampling{2,1})
    measPara.aniso(2) = length(measPara.oversampling{1,1});
end

% phase anistropy
if(isempty(measPara.ftlenPhase) || measPara.ftlenPhase == 0 || measPara.ftlenPhase ~= round(measPara.baseResolution * FOVphasePercentage * (1+phaseOversampling)))
    measPara.ftlenPhase = round(measPara.baseResolution * FOVphasePercentage * (1+phaseOversampling)); 
    if(measPara.ftlenPhase > measPara.baseResolution && FOVphasePercentage < 1) % can occur that y>x
        measPara.ftlenPhase = measPara.baseResolution;
    end
end
if(drecksMDH.Geo.PhaseRes(1) ~= 1)
    measPara.aniso(1) = measPara.ftlenPhase;
end
% phase oversampling
if(isempty(measPara.imageHeight) || measPara.imageHeight == 0)
    diffSize = abs(measPara.ftlenPhase - round(FOVphasePercentage * measPara.baseResolution));
    measPara.imageHeight = measPara.ftlenPhase;
else
    diffSize = abs(measPara.imageHeight - round(FOVphasePercentage * measPara.baseResolution));
end
if(mod(diffSize,2) == 0)
    oversamplimits = [round(diffSize/2)+1, round(diffSize/2)];
else
    oversamplimits = [round(diffSize/2)+1, 0];
    oversamplimits(2) = abs(diffSize-(oversamplimits(1)-1));
end
measPara.oversampling{1,2} = oversamplimits(1):(measPara.aniso(1)-oversamplimits(2));


if(strcmp(measPara.dimension,'3D') || strcmp(measPara.dimension,'4D') ||  strcmp(measPara.dimension,'5D'))
    % slice anisotropy    
    if(isempty(measPara.ftlenPar) || measPara.ftlenPar == 0 || measPara.ftlenPar ~= round(drecksMDH.Geo.ImagesPerSlab * (1+sliceOversampling)))
        measPara.ftlenPar = round(drecksMDH.Geo.ImagesPerSlab * (1+sliceOversampling));  
    end
    if(drecksMDH.Geo.PhaseRes(2) ~= 1)
        measPara.aniso(3) = measPara.ftlenPar;
    end
    % slice oversampling
    diffSize = abs(measPara.ftlenPar - drecksMDH.Geo.ImagesPerSlab);
    if(mod(diffSize,2) == 0)
        oversamplimits = [round(diffSize/2)+1, round(diffSize/2)];
    else
        oversamplimits = [round(diffSize/2)+1, 0];
        oversamplimits(2) = abs(diffSize-(oversamplimits(1)-1));
    end
    measPara.oversampling{1,3} = oversamplimits(1):(measPara.aniso(3)-oversamplimits(2));
%     measPara.imageDepth = measPara.imagesPerSlab;  %% ATTENTION !!!!!!!!!!! depth does not support time yet
else
    measPara.aniso(3) = measPara.LCall(1);
    measPara.oversampling{1,3} = 1:measPara.LCall(1);
%     measPara.imageDepth = measPara.LCall(1);     %% ATTENTION !!!!!!!!!!! depth does not support time yet
end


end

