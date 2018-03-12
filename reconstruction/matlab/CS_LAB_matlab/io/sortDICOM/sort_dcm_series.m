function  sort_dcm_series(root_dir)
% Sorts dcm files into folders by copying them


curr_dir=pwd;

% addpath(genpath(fileparts(fileparts(curr_dir))));

if(nargin < 1)
    error('Give a DICOM directory');
end

chdir(root_dir);
is_dicom_dir=0;
if exist('dicom','dir')
    dicomstr = 'dicom';
elseif(exist('DICOM','dir'))
    dicomstr = 'DICOM';
else
    dicomstr = '';
end
if(~exist([root_dir,filesep,dicomstr,'_ALL'],'dir'))
    copyfile(dicomstr, [dicomstr,'_ALL']);
end
if(~isempty(dicomstr))
    chdir(dicomstr);
    is_dicom_dir=1;
    %mkdir('../dicom_sorted');
    mkdir('../../dicom_sorted');
end
if(~exist('sort_dcm_series','file')) % sort_dcm_series not on path
    dicomMfilePath = fileparts(mfilename('fullpath'));
    addpath(genpath(dicomMfilePath));
end
files=get_file_list('*.IMA');
%files=get_file_list('*');

for f=1:length(files)
    dcm=dicominfo(files{f});
    dcm_path=[ regexprep(dcm.SeriesDescription,' ','_') '_' sprintf('%04d',dcm.SeriesNumber)];
    dcm_path=regexprep(dcm_path,'[*<>/&%$=+]','');
    if is_dicom_dir
        %dcm_path=[ '../dicom_sorted/' dcm_path];
        dcm_path=[ '../../dicom_sorted/' dcm_path];
    end
    if ~exist(dcm_path,'dir')
        mkdir(dcm_path);
        disp(dcm_path);
    end %if
    if exist(dcm_path,'dir')
        mv(files{f},dcm_path);
    end
 
end %f

chdir(root_dir);
mv('../dicom_sorted','.'); 
if(length(dir('DICOM')) == 2); rmdir DICOM; end
%status = rmdir('folderName')

chdir(curr_dir);

end
