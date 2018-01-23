function measPara = fMeasExtractMeasPara( drecksMDH, measPara, iLC, iLCPositions )
% extract meas para file
% input:
% drecksMDH     drecksMDH matrix or string containing path to MDH file
% measPara      struct for measurement parameters
% iLC           MDH loop counters
% iLCPositions  FOV and table positions
%
% output:
% measPara      struct with measurement parameters
%
% Copyright 2013/2014 Thomas Kuestner, University of Tuebingen, Germany
% thomas.kuestner@med.uni-tuebingen.de

if(ischar(drecksMDH))
    [drecksMDH, iLCPositions] = fMeasReadDrecksMDH(drecksMDH);
end

tmp = {'2D', '3D'};
measPara.dimensionSeq    = tmp{drecksMDH(1,1)+1};
if(drecksMDH(1,1) == 0) % 2D
    measPara.slices     = drecksMDH(1,2);
    measPara.partitions = 0;
else % 3D
    measPara.slices     = 0;
    measPara.partitions = drecksMDH(1,2);
end
if(drecksMDH(1,3) == 1) % MSM: sequential
    measPara.multisliceMode = 'sequential'; 
elseif(drecksMDH(1,3) == 2)
    measPara.multisliceMode = 'interleaved';
else
    measPara.multisliceMode = ''; % for 3D
end
% 1st MDH
measPara.multisliceModeInterleavingValue = drecksMDH(1,4);
measPara.sliceThickness     = drecksMDH(1,5)/100;
measPara.sliceDistance      = drecksMDH(1,6)/100;
measPara.lines              = drecksMDH(1,7);
measPara.baseResolution     = drecksMDH(1,8);
measPara.position           = ([double(drecksMDH(1,9)); double(drecksMDH(1,10)); double(drecksMDH(1,11))] - 20000)./10;
measPara.FOVread            = double(drecksMDH(1,12));
measPara.FOVphase           = double(drecksMDH(1,13));
measPara.phaseResolution    = drecksMDH(1,14)/100;
measPara.sliceResolution    = drecksMDH(1,16)/100;
measPara.imagesPerSlab      = drecksMDH(1,17);
measPara.TR                 = drecksMDH(1,18)/1000;
measPara.TE                 = drecksMDH(1,19)/1000;
measPara.flipAngle          = drecksMDH(1,20);
measPara.bandwidthPerPixel  = drecksMDH(1,21);
measPara.RFduration         = drecksMDH(1,22);
measPara.RFTimeBandwidth    = drecksMDH(1,23)/10;
measPara.CSAcceleration     = drecksMDH(1,24)/10;
measPara.CSFullySampled     = drecksMDH(1,25)/10;
% 2nd MDH
measPara.CSAcquiredLines    = drecksMDH(2,1);
measPara.CSFixedLines       = drecksMDH(2,2);
tmp = {'Poisson', 'Random', 'Gaussian'};
if(drecksMDH(2,3) && drecksMDH(2,3)-5 > 0 && drecksMDH(2,3)-5 < 4)
    measPara.CSSamplingType     = tmp{drecksMDH(2,3)-5};
else
    measPara.CSSamplingType     = 'Poisson';
end
tmp = {'None', 'Point', 'Block', 'Ellipse', 'Ring'};
if(drecksMDH(2,4) && drecksMDH(2,4) > 0 && drecksMDH(2,4) < 6)
    measPara.CSVDmap            = tmp{drecksMDH(2,4)};
else
    measPara.CSVDmap            = 'Ellipse';
end
measPara.CSaddPara          = [drecksMDH(2,5)/100, drecksMDH(2,6)];
measPara.averages           = drecksMDH(2,7);
measPara.concatenations     = drecksMDH(2,8);
measPara.sliceOversampling  = round(drecksMDH(2,9))/1000;    % round to 1st decimal place
measPara.phaseOversampling  = round(drecksMDH(2,10)/10)/100; % round to integer percentage values
measPara.readoutOversampling = drecksMDH(2,11);
tmp = {'y', 'z', 'off'};
measPara.ESPReSSoDirection  = tmp{drecksMDH(2,12)};
tmp = [81, 90, 99, 108, 117; 1, 0.5, 0.625, 0.75, 0.875];
measPara.ESPReSSo           = tmp(2,tmp(1,:) == drecksMDH(2,13));
if(strcmp(measPara.dimensionSeq,'2D') && strcmp(measPara.ESPReSSoDirection,'z')) % accidentally slice espresso was set in prior 3D scan
    measPara.ESPReSSoDirection = 'off';
    measPara.ESPReSSo       = 1;
end
measPara.imageHeight        = double(drecksMDH(2,14));
measPara.imageWidth         = double(drecksMDH(2,16));
measPara.FOVphasePercentage = measPara.imageHeight/measPara.imageWidth;
if(isnan(measPara.FOVphasePercentage))
    measPara.FOVphasePercentage = measPara.FOVphase/measPara.FOVread;
end
measPara.kspacesamplesPar   = drecksMDH(2,17);
measPara.ftlenPar           = drecksMDH(2,18);
measPara.kspacesamplesPhase = drecksMDH(2,19);
measPara.ftlenPhase         = drecksMDH(2,20);
measPara.echoLine           = drecksMDH(2,21);
measPara.echoPartition      = drecksMDH(2,22);
if(drecksMDH(2,23) == 42)
    measPara.sequenceName = 'CS_EPI';
elseif(drecksMDH(2,23) == 21)
    measPara.sequenceName = 'CS_Trufi';
elseif(drecksMDH(2,23) == 13)
    measPara.sequenceName = 'CS_Retro';
elseif(drecksMDH(2,23) == 1)
    measPara.sequenceName = 'CS_FLASH';
else
    measPara.sequenceName = 'CS_FLASH';
end
if(isempty(measPara.echoLine) || measPara.echoLine == 0)
    measPara.echoLine = round(measPara.lines/2);
end
center = zeros(1,3);
if(strcmp(measPara.dimensionSeq,'2D'))
    measPara.echoPartition = 0;
    if(mod(measPara.slices,2) == 0)
        measPara.echoSlice = measPara.slices/2 + 0.5;
    else
        measPara.echoSlice = ceil(measPara.slices/2);
    end
    center(3) = double(measPara.echoSlice);
    distSpacing = double(measPara.sliceDistance);
else % 3D
    measPara.echoSlice = 0;
    if(isempty(measPara.echoPartition) || measPara.echoPartition == 0)
        measPara.echoPartition = round(measPara.partitions/2);
    end 
    if(mod(measPara.echoPartition,2) == 0)
        center(3) = double(measPara.partitions/2) + 0.5;
    else
        center(3) = ceil(double(measPara.partitions/2));
    end
    distSpacing = double(measPara.sliceThickness);
end

% float        
measPara.normVector         = drecksMDH(:,15);
if(all(measPara.normVector == [1;0;0]))
    measPara.orientation    = 'sagittal';
elseif(all(measPara.normVector == [0;1;0]))
    measPara.orientation    = 'coronal';
elseif(all(measPara.normVector == [0;0;1]))
    measPara.orientation    = 'transversal';
else
    measPara.orientation    = 'rotated';
end

% correct calculated parameters (could be incorrect due to CS,
% ESPReSSo, EPI ref scans or rounding errors)
if(~isfield(measPara,'dim'))
    measPara.dim = ones(1,5);
    measPara.LCall = [0 1 1 1];
end
measPara.dim(1) = measPara.lines;
measPara.dim(2) = measPara.readoutOversampling * measPara.baseResolution; % due to CS_Trufi prescans
measPara.dim(3) = measPara.partitions;
if(strcmp(measPara.sequenceName, 'CS_Trufi') && ~isempty(iLC)) % number of kernels (SetLC) = cardiac phases = time domain
    measPara.dim(4) = double(max(iLC(:,10)) - min(iLC(:,10)) + 1); % SetLC
end
measPara.LCall(2) = measPara.averages;
if(drecksMDH(2,24) > 0)
    measPara.dim(4) = drecksMDH(2,24);
end
if(measPara.dim(3) <= 1)
    measPara.dimension = '2D';
else
    if(measPara.dim(4) > 1)
        measPara.dimension = '4D';
    else
        measPara.dimension = '3D';
    end    
end
if(drecksMDH(3,24) > 0)
    measPara.dim(5) = drecksMDH(3,24);
end
if(strcmp(measPara.sequenceName,'CS_Retro') && drecksMDH(2,25) > 0)
    measPara.dimension = '5D';
    measPara.dim(6) = drecksMDH(2,25);
else
    measPara.dim(6) = 1;
end

% initialize anisotropic correction (y-x-z/slice-t)
measPara.aniso = [measPara.dim(1), measPara.dim(2), measPara.dim(3), measPara.dim(4), measPara.dim(6)];
% oversampling correction exactly possible => no case switch necessary later
if(~isfield(measPara,'oversampling'))
    measPara.oversampling = cell(2,3);
end
measPara.oversampling{2,2} = true;
measPara.oversampling{2,3} = true;

% readout oversampling
readoutDiff = (measPara.readoutOversampling - 1) * measPara.baseResolution;
measPara.oversampling{1,1} = readoutDiff/2:measPara.dim(2)-readoutDiff/2-1;
if(isempty(measPara.imageWidth) || measPara.imageWidth == 0)
    measPara.imageWidth = length(measPara.oversampling{1,1});
end

% correct anisotropy if oversampling correction was done during recon
if(measPara.oversampling{2,1})
    measPara.aniso(2) = length(measPara.oversampling{1,1});
end

% phase anistropy
if(isempty(measPara.ftlenPhase) || measPara.ftlenPhase == 0 || measPara.ftlenPhase ~= round(measPara.baseResolution * measPara.FOVphasePercentage * (1+measPara.phaseOversampling)))
    measPara.ftlenPhase = round(measPara.baseResolution * measPara.FOVphasePercentage * (1+measPara.phaseOversampling)); 
    if(measPara.ftlenPhase > measPara.baseResolution)
        measPara.ftlenPhase = measPara.baseResolution;
    end
end
if(measPara.phaseResolution ~= 1)
    measPara.aniso(1) = measPara.ftlenPhase;
end
% phase oversampling
if(isempty(measPara.imageHeight) || measPara.imageHeight == 0)
    diffSize = abs(measPara.ftlenPhase - round(measPara.FOVphasePercentage * measPara.baseResolution));
    measPara.imageHeight = measPara.ftlenPhase;
else
    diffSize = abs(measPara.imageHeight - round(measPara.FOVphasePercentage * measPara.baseResolution));
end
if(mod(diffSize,2) == 0)
    oversamplimits = [round(diffSize/2)+1, round(diffSize/2)];
else
    oversamplimits = [round(diffSize/2)+1, 0];
    oversamplimits(2) = abs(diffSize-(oversamplimits(1)-1));
end
measPara.oversampling{1,2} = oversamplimits(1):(measPara.aniso(1)-oversamplimits(2));


if(strcmp(measPara.dimension,'3D') || strcmp(measPara.dimension,'4D'))
    % slice anisotropy    
    if(isempty(measPara.ftlenPar) || measPara.ftlenPar == 0 || measPara.ftlenPar ~= round(measPara.imagesPerSlab * (1+measPara.sliceOversampling)))
        measPara.ftlenPar = round(measPara.imagesPerSlab * (1+measPara.sliceOversampling));  
    end
    if(measPara.sliceResolution ~= 1)
        measPara.aniso(3) = measPara.ftlenPar;
    end
    % slice oversampling
    diffSize = abs(measPara.ftlenPar - measPara.imagesPerSlab);
    if(mod(diffSize,2) == 0)
        oversamplimits = [round(diffSize/2)+1, round(diffSize/2)];
    else
        oversamplimits = [round(diffSize/2)+1, 0];
        oversamplimits(2) = abs(diffSize-(oversamplimits(1)-1));
    end
    measPara.oversampling{1,3} = oversamplimits(1):(measPara.aniso(3)-oversamplimits(2));
    measPara.imageDepth = double(measPara.imagesPerSlab);  %% ATTENTION !!!!!!!!!!! depth does not support time yet
else
    measPara.aniso(3) = measPara.LCall(1);
    measPara.oversampling{1,3} = 1:measPara.LCall(1);
    measPara.imageDepth = double(measPara.LCall(1));     %% ATTENTION !!!!!!!!!!! depth does not support time yet
end

% calculate voxel positions [mm] in global scanner coordinate system
%           A
%           |      H  
%           |     ^ z
%           |    /
%           |   /
%           |  /
%           | /
%           |/
% R --------|-------------> L  x
%          /|
%         / |
%        F  |
%           v  y
%           P
%
% WATCH OUT: Unknown offset between DICOM and MDH position! MDH
% position was verified and can be trusted!

% get the pixel spacing and central k-space sample (y-x-z/slice)
measPara.pixelSpacing = [measPara.FOVphase/measPara.imageHeight, measPara.FOVread/measPara.imageWidth, distSpacing];
center(1:2) = [measPara.imageHeight/2, measPara.imageWidth/2];   

switch measPara.orientation
    case 'sagittal'
        measPara.xPos = measPara.position(1) + ([1:measPara.imageDepth].' - center(3)) .* measPara.pixelSpacing(3);
        measPara.yPos = measPara.position(2) + ([1:measPara.imageHeight].' - center(1)) .* measPara.pixelSpacing(1);
        measPara.zPos = measPara.position(3) + ([1:measPara.imageWidth].' - center(2)) .* measPara.pixelSpacing(2);
    case 'coronal'
        measPara.xPos = measPara.position(1) + ([1:measPara.imageWidth].' - center(2)) .* measPara.pixelSpacing(2);
        measPara.yPos = measPara.position(2) + ([1:measPara.imageDepth].' - center(3)) .* measPara.pixelSpacing(3);
        measPara.zPos = measPara.position(3) + ([1:measPara.imageHeight].' - center(1)) .* measPara.pixelSpacing(1);
    case 'transversal'
        measPara.xPos = measPara.position(1) + ([1:measPara.imageWidth].' - center(2)) .* measPara.pixelSpacing(2);
        measPara.yPos = measPara.position(2) + ([1:measPara.imageHeight].' - center(1)) .* measPara.pixelSpacing(1);
        measPara.zPos = measPara.position(3) + ([1:measPara.imageDepth].' - center(3)) .* measPara.pixelSpacing(3);
    otherwise 
        %
end


if(exist('iLCPositions','var') && ~isempty(iLCPositions))
    % additional parameters
    measPara.sliceDataPosition = iLCPositions(:,1:3);
    measPara.PTAB              = -((iLCPositions(:,4) +5)/10);  
    measPara.rotMatrix         = quat2dcm(iLCPositions(1,5:8));
    measPara.rotMatrix(:,2:3)  = -measPara.rotMatrix(:,2:3); % due to negative orientation of y and z axis (compared to standard cartesian system)
    [yaw, pitch, roll]         = dcm2angle(measPara.rotMatrix);
    measPara.rotAngles         = [yaw, pitch, roll];
end


end

