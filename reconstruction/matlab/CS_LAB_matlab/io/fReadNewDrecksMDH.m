function SInfo = fReadNewDrecksMDH(dMDH)
% read and process drecksMDH v2.0
%
% (c) Christian Wuerslin, Thomas Kuestner, 2013
% ---------------------------------------------------------------------

iPos = 1;

while iPos <= length(dMDH)
    iTag = uint16(dMDH(iPos));
    iPos = iPos + 1;
    lTag = bitget(iTag, 16);
    while ~lTag
        if iPos > length(dMDH), break, end
        iTag = uint16(dMDH(iPos));
        iPos = iPos + 1;
        lTag = bitget(iTag, 16);
    end
    if iPos > length(dMDH), iPos = length(dMDH); end
    iLength = bitand(iTag, 7);
    iPower = bitand(iTag, 56)/8;
    iID = bitand(iTag, 1984)/64;
    iGroup = bitand(iTag, 30720)/2048;
    dData = dMDH(iPos:(iPos + iLength));
    dData = dData./(10^double(iPower));
    iPos = iPos + iLength + 1;
    
    switch iGroup
        case 01
            sGroup = 'Seq';
            switch iID
                case 01
                    sID = 'Sequence';
                    switch dData
                        case 1,  dData = 'CS_FLASH';
                        case 13, dData = 'CS_Retro';
                        case 21, dData = 'CS_Trufi';
                        case 42, dData = 'CS_EPI';
                        otherwise
                    end
                case 02
                    sID = 'Version'; 
                    dData = sprintf('%d.%d', dData(:)');
                case 03
                    sID = 'Build';
                    dData = sprintf('%d.%d.%d %02d:%02d:%02d', dData(:)');
                case 04, sID = 'MSM';
                case 05, sID = 'Bandwidth';
                case 06, sID = 'Is3D';
                otherwise, sID = sprintf('ID%02d', iID);
            end
        case 02
            sGroup = 'LC';
            switch iID
                case 01, sID = 'Measurements';
                case 02, sID = 'Averages';
                case 03, sID = 'Concats';
                otherwise, sID = sprintf('ID%02d', iID);
            end
        case 03
            sGroup = 'Geo';
            switch iID
                case 01, sID = 'Shift'; dData = dData - 2000;
                case 02, sID = 'MatrixSize';
                case 03, sID = 'NormSag';
                case 04, sID = 'NormCor';
                case 05, sID = 'NormTra';
                case 06, sID = 'FOV';
                case 07, sID = 'PhaseRes';
                case 08, sID = 'ImagesPerSlab';
                case 09, sID = 'OverSampling';
                case 10, sID = 'Spacing';
                case 11, sID = 'ImageSize';
                case 12, sID = 'EchoPosition';
                case 13, sID = 'FFTLength';
                otherwise, sID = sprintf('ID%02d', iID);
            end
        case 04
            sGroup = 'Contrast';
            switch iID
                case 01, sID = 'TR';
                case 02, sID = 'TE';
                case 03, sID = 'FlipAngle';
                otherwise, sID = sprintf('ID%02d', iID);
            end
        case 05
            sGroup = 'Accel';
            switch iID
                case 01, sID = 'EspressoDir';
                case 02, sID = 'EspressoFactor';
                otherwise, sID = sprintf('ID%02d', iID);
            end
        case 14
            sGroup = 'Wip';
            if(strcmp(SInfo.Seq.Sequence,'CS_Retro') && str2double(SInfo.Seq.Version) < 2.1) % old version, backward compatibility
                switch iID
                    case 01, sID = 'RFDuration';
                    case 02, sID = 'TimeBandwidth';
                    case 03, sID = 'NavPeriod';
                    case 04, sID = 'NavRes';
                    case 05, sID = 'SamplingFactor';
                    case 06, sID = 'WFactor';
                    case 07, sID = 'Phases';
                    case 08, sID = 'TotalScanTime';
                    case 09, sID = 'MotionModelTime';    
                    case 10, sID = 'NMotionModel';                                    
                    case 11, sID = 'RespPeriod'; 
                    case 12, sID = 'ECGPeriod';                    
                    otherwise, sID = sprintf('ID%02d', iID);
                end
            elseif(strcmp(SInfo.Seq.Sequence,'CS_Retro') && str2double(SInfo.Seq.Version) == 2.1) % old version, backward compatibility
                switch iID
                    case 01, sID = 'RFDuration';
                    case 02, sID = 'TimeBandwidth';
                    case 03, sID = 'NavPeriod';
                    case 04, sID = 'NavRes';
                    case 05, sID = 'SamplingFactor';
                    case 06, sID = 'WFactor';
                    case 07, sID = 'Phases';
                    case 08, sID = 'TotalScanTime';
                    case 09, sID = 'UNKNOWN';    
                    case 10, sID = 'Mask';                                    
                    case 11, sID = 'MotionModelTime';
                    case 12, sID = 'NMotionModel';          
                    case 13, sID = 'RespPeriod';
                    case 14, sID = 'ECGPeriod';
                    otherwise, sID = sprintf('ID%02d', iID); 
                end
            elseif(strcmp(SInfo.Seq.Sequence,'CS_FLASH') && str2double(SInfo.Seq.Version) < 3)
                switch iID
                    case 01, sID = 'RFDuration';
                    case 02, sID = 'TimeBandwidth';
                    case 03, sID = 'NavPeriod';
                    case 04, sID = 'NavRes';
                    case 05, sID = 'SamplingFactor';
                    case 06, sID = 'WFactor';
                    case 07, sID = 'Phases';
                    case 08, sID = 'TotalScanTime';
                    case 09, sID = 'FullySampled';
                    case 10
                        sID = 'Mask';
                        dData = sprintf('%02d,%02d', dData(:)');
                    otherwise, sID = sprintf('ID%02d', iID);
                end
            else % new (common) version
                switch iID
                    case 01, sID = 'RFDuration';
                    case 02, sID = 'TimeBandwidth';
                    case 03, sID = 'NavPeriod';
                    case 04, sID = 'NavRes';
                    case 05, sID = 'SamplingFactor';
                    case 06, sID = 'WFactor';
                    case 07, sID = 'Phases';
                    case 08, sID = 'TotalScanTime';
                    case 09, sID = 'MotionModelTime';
                    case 10, sID = 'NMotionModel';
                    case 11, sID = 'RespPeriod';
                    case 12, sID = 'ECGPeriod';
                    case 13
                        sID = 'Mask';
                        dData = sprintf('%02d,%02d', dData(:)');
                    case 14, sID = 'FullySampled';
                    case 15, sID = 'ShortTrajAccel';
                    case 16, sID = 'LongTrajAccel';
                    case 17, sID = 'NSpiralTurn';
                    otherwise, sID = sprintf('ID%02d', iID);
                end
            end
        case 15
            sGroup = 'Misc';
            sID = sprintf('ID%02d', iID);
        otherwise
            sGroup = sprintf('Group%02d', iGroup);
            sID = sprintf('ID%02d', iID);
    end
    eval(sprintf('SInfo.%s.%s = dData;', sGroup, sID));
end