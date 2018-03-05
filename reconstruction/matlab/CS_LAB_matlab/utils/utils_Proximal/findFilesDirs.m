function [ loadFiles, loadDirs ] = findFilesDirs( imagePath)
% recursive search
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

[filesInPath, dirsInPath] = rdir(imagePath);
% filter files (just use *.mat and *.dat)
allFiles = cell(1,length(filesInPath));
allDirs = cell(1,length(dirsInPath));
foundFile = false(1,length(filesInPath));
if ~isempty(dirsInPath)
    nDirsInPath = length(dirsInPath);
    for k = 1:nDirsInPath
        allDirs{k} = dirsInPath.fullpath;
    end;
else
    nDirsInPath = 1;
end;
foundDir = false(1,length(nDirsInPath));

for i=1:length(filesInPath)
    allFiles{i} = filesInPath(i).fullpath;
        tmp = regexp(allFiles{i},'\w+x\d*\w+\.[md]at', 'once');
        tmp2 = regexp(allFiles{i},'\w+_rec\w+', 'once');
        if (~isempty(tmp) && isempty(tmp2))
            foundFile(i) = true;
        end;
end;
for k=1:nDirsInPath; % for DICOM only check for dirs, since fReadDICOM searches a whole folder -> if folder contains *.IMA-file, set found(k) true
    for i=1:length(filesInPath)
        if ~isempty(allDirs)
            tmp = regexp(allFiles{i},'\w+.IMA', 'once');
            if(~isempty(tmp))
                dirString = regexptranslate('escape', allDirs{k});
                tmp_dir = regexp(allFiles{i},[dirString '\w*'],'once');
                if(~isempty(tmp_dir))
                    foundDir(k) = true;
                end;
            end;
        else
            allDirs{k} = imagePath;
            allFiles{i} = filesInPath(i).fullpath;
            tmp = regexp(allFiles{i},'\w+.IMA', 'once');
            if(~isempty(tmp))
                dirString = regexptranslate('escape', allDirs{k});
                regexp(allFiles{i},[dirString '\w*']);
                foundDir(k) = true;
            end;
        end;
    end;
end

loadFiles = allFiles(foundFile);
loadDirs = allDirs(foundDir);

% remove same file names (only needed for mat/dat)
helper = cellfun(@(x) x(1:end-4), loadFiles, 'UniformOutput', false); % removes file extension in string
[~, idxloadFiles, ~] = unique(helper); % sorts out double entries
loadFiles = loadFiles(idxloadFiles);

end
