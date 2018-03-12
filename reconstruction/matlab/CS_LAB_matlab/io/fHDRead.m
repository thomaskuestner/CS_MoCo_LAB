function xData = fHDRead(sFilename)
% read in uncompacted mat file
%
% (c) Christian Wuerslin, Thomas Kuestner
% ---------------------------------------------------------------------

fid = fopen(sFilename, 'r');
sCtlStr = char(fread(fid, 5, 'char'));
if ~strcmp(sCtlStr', 'matHD'), error('Not a matHD file.'); end

lComplex = false;
iFormat = fread(fid, 1, 'uint8');
if iFormat >= 16
    lComplex = true;
    iFormat = iFormat - 16;
end
switch(iFormat)
    case  8, sClass = 'uint8';
    case  9, sClass = 'uint16';
    case 10, sClass = 'uint32';
    case 11, sClass = 'uint64';
    case 12, sClass = 'int8';
    case 13, sClass = 'int16';
    case 14, sClass = 'int32';
    case 15, sClass = 'int64';
    case  2, sClass = 'single';
    case  3, sClass = 'double';
    case  1, sClass = 'logical';
    otherwise
        error('Data type not supported.');
end

iNDims = fread(fid, 1, 'uint16');
iSize  = uint16(fread(fid, iNDims, 'uint16'));

xData = fread(fid, prod(iSize), ['*', sClass]);
if lComplex
    xData = complex(xData, fread(fid, prod(iSize), ['*', sClass]));
end
fclose(fid);

if iFormat == 1, xData = logical(xData); end

xData = reshape(xData, iSize');