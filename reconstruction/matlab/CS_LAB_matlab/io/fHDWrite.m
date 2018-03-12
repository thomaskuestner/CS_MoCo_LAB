function fHDWrite(sFilename, xData)
% write uncompacted mat file
%
% (c) Christian Wuerslin, Thomas Kuestner
% ---------------------------------------------------------------------

switch(class(xData))
    case 'uint8' , iFormat =  8;
    case 'uint16', iFormat =  9;
    case 'uint32', iFormat = 10;
    case 'uint64', iFormat = 11;
    case 'int8'  , iFormat = 12;
    case 'int16' , iFormat = 13;
    case 'int32' , iFormat = 14;
    case 'int64' , iFormat = 15;
    case 'single', iFormat =  2;
    case 'double', iFormat =  3;
    case 'logical'
        iFormat = 1;
        xData = uint8(xData);
    
    otherwise
        error('Data type not supported.');
end

if ~isreal(xData), iFormat = iFormat + 16; end

fid = fopen(sFilename, 'w');
fwrite(fid, 'matHD', 'char');
fwrite(fid, uint8(iFormat), 'uint8');
fwrite(fid, uint16(ndims(xData)), 'uint16');
fwrite(fid, uint16(size(xData)), 'uint16');
if isreal(xData)
    fwrite(fid, xData, class(xData));
else
    fwrite(fid, real(xData), class(xData));
    fwrite(fid, imag(xData), class(xData));
end
fclose(fid);