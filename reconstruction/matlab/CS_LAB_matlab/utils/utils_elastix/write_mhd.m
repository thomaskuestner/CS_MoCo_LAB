function f = write_mhd(filename, image, varargin)
% Write mhd
% arguments:
% write_mhd(filename,image)
% write_mhd(filename,image, param, value,...)
% params and default value:
%   'NDims'                     =   3
%   'BinaryData'                =   true
%   'BinaryDataByteOrderMSB'    =   false
%   'CompressedData'            =   false
%   'TransformMatrix'           =   1 0 0 0 1 0 0 0 1 (coherent with NDims)
%   'CenterOfRotation'          =   0 0 0 (coherent with NDims)
%   'AnatomicalOrientation'     =   RAI
%   'ElementNumberOfChannels'   =   1
%   'ElementType'               =   MET_FLOAT
%
% modifications by Thomas Kuestner, 2015
% ---------------------------------------------------------------------


% revise this
element_types=struct('uchar','MET_UCHAR','double','MET_DOUBLE','single','MET_FLOAT','logical','MET_CHAR','int8','MET_UCHAR','uint8','MET_UCHAR','int16','MET_SHORT','uint16','MET_USHORT','int32','MET_INT','uint32','MET_UINT');
% element_types=struct('double','MET_DOUBLE','int8','MET_CHAR','uint8','MET_UCHAR','int16','MET_SHORT','uint16','MET_USHORT','int32','MET_INT','uint32','MET_UINT');

%Look for arguments
NDims                       = numel(image.size);
BinaryData                  = 'True';
BinaryDataByteOrderMSB      = 'False';
CompressedData              = 'False';

TransformMatrix  = num2str(reshape(image.orientation,1,NDims*NDims)); %'1 0 0 0 1 0 0 0 1';% (coherent with NDims)
CenterOfRotation = num2str(zeros(1,NDims)) ;

AnatomicalOrientation =   'RAI';
if (isa(image,'VectorImageType') || isfield(image,'datax'))
    ElementNumberOfChannels = 3;
else
    ElementNumberOfChannels = 1;
end


%ElementType =getfield(element_types,class(image.data));
EType ='double';


for i=1:size(varargin,2)
    switch (lower(varargin{i}))
        case 'ndims'
            NDims =  varargin{i+1};
        case 'binarydata'
            BinaryData  =   varargin{i+1};
        case 'binarydatabyteordermsb'
            BinaryDataByteOrderMSB =   varargin{i+1};
        case 'compresseddata'
            CompressedData =   varargin{i+1};
        case 'transformmatrix'
            TransformMatrix =   varargin{i+1};
        case 'centerofrotation'
            CenterOfRotation =   varargin{i+1};
        case 'anatomicalorientation'
            AnatomicalOrientation =   varargin{i+1};
        case 'elementnumberofchannels'
            ElementNumberOfChannels =   varargin{i+1};
        case 'elementtype'
            EType =  varargin{i+1};
    end
end
ElementType =  element_types.(EType );
% Extract path and filenames---
path_ = regexp(filename,'/','split');
path = [];
for i=2:(size(path_,2)-1)
     
    path = [path '/' path_{i}];
end

raw_ = regexp(filename,'\.','split');
raw = [];
for i=1:(size(raw_,2)-1)
    raw = [raw  raw_{i} '.'];
end
rawfile = [raw 'raw'];
rawfilename_ = regexp(rawfile,'/','split');
rawfilename = rawfilename_{size(rawfilename_,2)};
% write *.mhd file
fid=fopen(filename,'w','native');
fprintf(fid,'ObjectType = Image\n');
fprintf(fid,'NDims = %d\n',NDims);
fprintf(fid,'BinaryData = %s\n',BinaryData);
fprintf(fid,'BinaryDataByteOrderMSB = %s\n',BinaryDataByteOrderMSB);
fprintf(fid,'CompressedData = %s\n',CompressedData);
fprintf(fid,'TransformMatrix = %s\n',TransformMatrix);
fprintf(fid,'CenterOfRotation = %s\n',CenterOfRotation);
fprintf(fid,'AnatomicalOrientation = %s\n',AnatomicalOrientation);
fprintf(fid,'%s\n',['Offset = ' num2str(image.origin(:)') ]);
fprintf(fid,'%s\n',['ElementSpacing = ' num2str(image.spacing(:)') ]);
fprintf(fid,'%s\n',['DimSize = ' num2str(image.size(:)') ]);

fprintf(fid,'ElementNumberOfChannels = %d\n',ElementNumberOfChannels);
fprintf(fid,'ElementType = %s\n',ElementType);
fprintf(fid,'ElementDataFile = %s\n',rawfilename);

fclose(fid);

% write *.raw file

%dataToWrite=image.data;

%if (islogical(dataToWrite))
%     dataToWrite = cast(dataToWrite,'uint8');
% else

    
    
fid=fopen(rawfile,'w','native');
if (ElementNumberOfChannels ==1)
    %fwrite(fid,image.data,class(dataToWrite) );
    dataAll = cast(image.data,EType);    
elseif (ElementNumberOfChannels==2)
    %do nothing
elseif (ElementNumberOfChannels==3)
    %write point per point
    dataAll = cast(permute(cat(NDims+1,image.datax,  image.datay, image.dataz),...
        [NDims+1 1:NDims]),EType);
elseif (ElementNumberOfChannels==4)
    %TODO
end
fwrite(fid,dataAll,EType);
fclose(fid);

f=fid;

end