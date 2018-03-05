function [ samplingMask ] = readSamplingMaskFromFile(sFilename)
%This function reads out a sampling pattern from a .txt-file

if(nargin < 1)
    file = fopen('samplingPattern.txt','r');
else
    file = fopen(sFilename,'r');
end

endOfDoc = false;
i = 0;

while(~endOfDoc)
    actualLine = fgetl(file);
    
    if(actualLine == -1)
        endOfDoc = true;
        break;
    end
    i = i + 1;
    lines{i} = actualLine;
    
    for j=1:numel(lines{i})
        actualChar = lines{i}(j);

        if(actualChar == '-')
            samplingMask(i, j) = 0;
        elseif (actualChar == 'x' || actualChar == '*' || actualChar == 'O' || actualChar == '0')
            samplingMask(i, j) = 1;
        else
            error('not working, because there is an unknown charakter');
        end
    end
end

if(exist('drawMask','file'))    
    drawMask(samplingMask);
end

fclose(file);

end

