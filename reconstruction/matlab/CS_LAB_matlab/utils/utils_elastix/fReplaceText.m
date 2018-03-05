function fReplaceText(fFilename, sText, sWhatReplaced)
% find and replace text in txt file
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

    if(nargin < 3)
        sWhatReplaced = 'InitialTransformParametersFileName';
    end
    fid = fopen(fFilename,'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
        if(ischar(A{i}) && ~isempty(regexp(A{i},[sprintf('^(%s',sWhatReplaced),'\w*'],'once')))
            A{i} = sprintf('(%s "%s")',sWhatReplaced,sText);
        end
    end
    fclose(fid);
    fid = fopen(fFilename, 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break;
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);
end