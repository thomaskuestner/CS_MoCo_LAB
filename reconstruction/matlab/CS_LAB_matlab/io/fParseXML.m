function SDrecksMDH  = fParseXML( sFilename, measOffset )
% parse XML file to drecksMDH 2.0

%sXML = xml2struct(sFilenameXML);

fidmeas = fopen(sFilename, 'r', 'ieee-le');
fseek(fidmeas,measOffset{2},'bof');
iHeaderLength = fread(fidmeas, 1, 'uint32');
        
[hdr, rstraj] = read_twix_hdr(fidmeas);

% Seq
SDrecksMDH.Seq.Sequence         = hdr.Meas.ProtocolName;
SDrecksMDH.Seq.Version          = '4.x';
SDrecksMDH.Seq.Build            = hdr.Meas.SequenceDescription;
SDrecksMDH.Seq.Bandwidth        = [];
if(strcmp(hdr.Dicom.tMRAcquisitionType, '3D'))
    SDrecksMDH.Seq.Is3D             = true;
else
    SDrecksMDH.Seq.Is3D             = false;
end
SDrecksMDH.Seq.MSM              = hdr.Meas.MultiSliceMode;
% LC
SDrecksMDH.LC.Measurements      = hdr.Meas.NRepMeas;
SDrecksMDH.LC.Averages          = hdr.Meas.NAve;
SDrecksMDH.LC.Concats           = hdr.Meas.NSlcMeas;
% Geo
SDrecksMDH.Geo.Shift            = [hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dSag,hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor,hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dTra];
SDrecksMDH.Geo.MatrixSize       = [hdr.MeasYaps.sKSpace.lBaseResolution, hdr.MeasYaps.sKSpace.lPhaseEncodingLines, hdr.MeasYaps.sKSpace.lPartitions];
if(isfield(hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dSag'))
    SDrecksMDH.Geo.NormSag      = hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dSag;
else
    SDrecksMDH.Geo.NormSag      = 0;
end
if(isfield(hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dCor'))
    SDrecksMDH.Geo.NormCor      = hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dCor;
else
    SDrecksMDH.Geo.NormCor      = 0;
end
if(isfield(hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dTra'))
    SDrecksMDH.Geo.NormTra      = hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dTra;
else
    SDrecksMDH.Geo.NormTra      = 0;
end
SDrecksMDH.Geo.FOV              = [hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV,hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV, hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness];
SDrecksMDH.Geo.PhaseRes         = hdr.MeasYaps.sKSpace.dPhaseResolution;
SDrecksMDH.Geo.ImagesPerSlab    = hdr.MeasYaps.sKSpace.lImagesPerSlab;
if(isempty(hdr.Dicom.flPhaseOS)), flPhaseOS = 0; else flPhaseOS = hdr.Dicom.flPhaseOS; end;
if(isempty(hdr.Dicom.flSliceOS)), flSliceOS = 0; else flSliceOS = hdr.Dicom.flSliceOS; end;
SDrecksMDH.Geo.Oversampling     = [hdr.Dicom.flReadoutOSFactor, flPhaseOS, flSliceOS];
SDrecksMDH.Geo.Spacing          = SDrecksMDH.Geo.FOV./SDrecksMDH.Geo.MatrixSize;
SDrecksMDH.Geo.ImageSize        = [hdr.Meas.NImageCols, hdr.Meas.NImageLins, hdr.Meas.NImagePar];
SDrecksMDH.Geo.EchoPosition     = ceil(SDrecksMDH.Geo.MatrixSize(2:3)./2);
SDrecksMDH.Geo.FFTLength        = [];
% Contrast
SDrecksMDH.Contrast.TR          = hdr.Meas.alTR(1)./1000;
SDrecksMDH.Contrast.TE          = hdr.Meas.alTE(1)./1000;
SDrecksMDH.Contrast.FlipAngle   = hdr.Meas.adFlipAngleDegree;
% Accel
switch hdr.MeasYaps.sWipMemBlock.alFree{3}
    case 1 % off
        SDrecksMDH.Accel.EspressoDir    = 'y';
        SDrecksMDH.Accel.EspressoFactor = 1;
    case 2 % y05
        SDrecksMDH.Accel.EspressoDir    = 'y';
        SDrecksMDH.Accel.EspressoFactor = 0.5;
    case 3 % y58
        SDrecksMDH.Accel.EspressoDir    = 'y';
        SDrecksMDH.Accel.EspressoFactor = 5/8;
    case 4 % y68
        SDrecksMDH.Accel.EspressoDir    = 'y';
        SDrecksMDH.Accel.EspressoFactor = 6/8;
    case 5 % y78
        SDrecksMDH.Accel.EspressoDir    = 'y';
        SDrecksMDH.Accel.EspressoFactor = 7/8;
    case 6 % z05
        SDrecksMDH.Accel.EspressoDir    = 'z';
        SDrecksMDH.Accel.EspressoFactor = 0.5;
    case 7 % z58
        SDrecksMDH.Accel.EspressoDir    = 'z';
        SDrecksMDH.Accel.EspressoFactor = 5/8;
    case 8 % z68
        SDrecksMDH.Accel.EspressoDir    = 'z';
        SDrecksMDH.Accel.EspressoFactor = 6/8;
    case 9 % z78
        SDrecksMDH.Accel.EspressoDir    = 'z';
        SDrecksMDH.Accel.EspressoFactor = 7/8;
end
% WIP
SDrecksMDH.Wip.RFDuration       = hdr.MeasYaps.sWipMemBlock.alFree{4};
SDrecksMDH.Wip.TimeBandwidth    = hdr.MeasYaps.sWipMemBlock.adFree{1};
SDrecksMDH.Wip.NavPeriod        = hdr.MeasYaps.sWipMemBlock.adFree{7};
SDrecksMDH.Wip.NavRes           = [hdr.MeasYaps.sWipMemBlock.alFree{5}, hdr.MeasYaps.sWipMemBlock.alFree{6}];
SDrecksMDH.Wip.SamplingFactor   = hdr.MeasYaps.sWipMemBlock.adFree{3};           
SDrecksMDH.Wip.TotalScanTime    = hdr.MeasYaps.sWipMemBlock.alFree{8};
SDrecksMDH.Wip.MotionModelTime  = hdr.MeasYaps.sWipMemBlock.alFree{2};
SDrecksMDH.Wip.NMotionModel     = hdr.MeasYaps.sWipMemBlock.alFree{10};
SDrecksMDH.Wip.RespPeriod       = hdr.MeasYaps.sWipMemBlock.adFree{8};
SDrecksMDH.Wip.ECGPeriod        = hdr.MeasYaps.sWipMemBlock.adFree{9};
SDrecksMDH.Wip.Mask             = hdr.MeasYaps.sWipMemBlock.alFree{12};
SDrecksMDH.Wip.Recon            = hdr.MeasYaps.sWipMemBlock.alFree{11};
SDrecksMDH.Wip.FullySampled     = hdr.MeasYaps.sWipMemBlock.adFree{2};
SDrecksMDH.Wip.ShortTrajAccel   = hdr.MeasYaps.sWipMemBlock.adFree{4};
SDrecksMDH.Wip.LongTrajAccel    = hdr.MeasYaps.sWipMemBlock.adFree{5};
SDrecksMDH.Wip.NSpiralTurn      = hdr.MeasYaps.sWipMemBlock.adFree{6};


end

