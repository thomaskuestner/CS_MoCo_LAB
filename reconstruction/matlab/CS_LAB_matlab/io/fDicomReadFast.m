function [iImg, iHdr] = fDicomReadFast(sFilename, iWidth, iHeight)
% read in DICOM files
%
% (c) Christian Wuerslin, Thomas Kuestner, 2011
% ---------------------------------------------------------------------

iWORDSIZE = 2;

iHeight = double(iHeight);
iWidth = double(iWidth);

SDir = dir(sFilename);
iBytesInFile = SDir(1).bytes;
iHeaderSize = iBytesInFile - iWORDSIZE*iWidth*iHeight;

fid = fopen(sFilename);
% fseek(fid, iHeaderSize, 'bof');
iHdr = fread(fid, iHeaderSize, '*uint8');

iImg = zeros(iHeight, iWidth, 'uint16');
iImg(:) = fread(fid, iWidth*iHeight, '*uint16');

fclose(fid);

iImg = transpose(iImg);