function SImg = fMCReconMain(dImg, SImg, sElastixPath, mcPath, mcParaPath, sElastixParamFile, bType, bCalcDF)
% main motion correction (MC) function
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

if(nargin < 7)
    bType = 0;
    bCalcDF = false;
end

if(nargin < 8)
    bCalcDF = false;
end

% bType: binary type 
% -> bInit      001 = 1
% -> bCompose   010 = 2
bInit = false;
bCompose = false; % bCompose implies bInit (since calling of elastix with "-t0" makes no sense for composed trafo) -> take same if-path as bInit
if(bType == 1)
    bInit = true;
elseif(bType == 2)
    bCompose = true;
    warning('fMCReconMain(): bInit should be true!');
elseif(bType == 3)
    bInit = true;
    bCompose = true;
end

% check for SURF features -> dependent on mcParaPath
[~,bSurf] = fileparts(mcParaPath);
if(strcmp(bSurf,'SURF')) 
    bSurf = true;
% elseif(strcmp(bSurf,'compose'))
%     bSurf = false; % currently not supported -> TODO: make variable for distinction
%     bCompose = true;
else % normal
    bSurf = false;
end

iNGATES = size(dImg,4);
if(~isempty(sElastixParamFile))
    for iMC = 1:iNGATES
%         dImg(:,:,:,iMC) = scaleImg(dImg(:,:,:,iMC),[0 2^16-1]);
%         SImg.data = uint16(dImg(:,:,:,iMC)); % TODO: check if scaling over complete image (including time) or (like here) just scaling over single resp states
%         write_mhd([mcPath,filesep,sprintf('Gate%02u.mhd', iMC)], SImg, 'ElementType', 'uint16');  
        
        SImg.data = dImg(:,:,:,iMC);
        write_mhd([mcPath,filesep,sprintf('Gate%02u.mhd', iMC)], SImg, 'ElementType', 'single');
    end
end

if(bCompose)
    sElastixParamFileInit = sElastixParamFile;
end

if(bSurf)
    iTotal = (iNGATES-1)*(size(dImg,3)+2);
else
    iTotal = (iNGATES-1)*2;
end
dispProgress('MC',0,iTotal);
for iMC = 2:iNGATES % TODO: evtl. parfor => einzelne Unterordner für jedes Gate erstellen, sonst werden Daten überschrieben
    % include SURF features
    if(~isempty(sElastixParamFile) && bSurf) % shortcut to just retrieve deformation field
        xFixed = cell(1,size(dImg,3));
        xMoved = xFixed;
        yFixed = xFixed;
        yMoved = xFixed;
        for iZ = 1:size(dImg,3)
            [fFixed,dFixed] = vl_sift(single(dImg(:,:,iZ,1)));
            [fMoved,dMoved] = vl_sift(single(dImg(:,:,iZ,iMC)));
            [matches, scores] = vl_ubcmatch(dFixed,dMoved,2);
            [~, perm] = sort(scores, 'descend');
            matches = matches(:, perm);
    %         scores  = scores(perm);
            xFixed{iZ} = fFixed(1,matches(1,:));
            xMoved{iZ} = fMoved(1,matches(2,:));
            yFixed{iZ} = fFixed(2,matches(1,:));
            yMoved{iZ} = fMoved(2,matches(2,:));
            dispProgress('MC',((iMC-2)*(size(dImg,3)+2)+iZ)/iTotal);
        end
        nPoints = sum(cellfun(@length, xFixed));
        fidFixed = fopen([mcPath,filesep,sprintf('PointSetGate%02uTo01_fixed.txt',iMC)],'w');
        fidMoved = fopen([mcPath,filesep,sprintf('PointSetGate%02uTo01_moved.txt',iMC)],'w');
        fprintf(fidFixed,'index\n%d\n',nPoints);
        fprintf(fidMoved,'index\n%d\n',nPoints);
        for iZ = 1:size(dImg,3)
            for iInner = 1:length(xFixed{iZ})
                if(iZ == size(dImg,3) && iInner == length(xFixed{iZ}))
                    fprintf(fidFixed,'%.4f %.4f %.4f', xFixed{iZ}(iInner)-1, yFixed{iZ}(iInner)-1, iZ-1);
                    fprintf(fidMoved,'%.4f %.4f %.4f', xMoved{iZ}(iInner)-1, yMoved{iZ}(iInner)-1, iZ-1);
                else
                    fprintf(fidFixed,'%.4f %.4f %.4f\n', xFixed{iZ}(iInner)-1, yFixed{iZ}(iInner)-1, iZ-1);
                    fprintf(fidMoved,'%.4f %.4f %.4f\n', xMoved{iZ}(iInner)-1, yMoved{iZ}(iInner)-1, iZ-1);
                end
            end
        end
        fclose(fidFixed);
        fclose(fidMoved);
    end
    if(bCompose)
        % adapt paramfile
        for iParams=1:length(sElastixParamFileInit)
            % forward
            copyfile([mcParaPath,filesep,sElastixParamFileInit{iParams}], [mcParaPath,filesep,sElastixParamFileInit{iParams}(1:end-6),sprintf('%02u.txt',iMC)],'f');
            sElastixParamFile{iParams} = [sElastixParamFile{iParams}(1:end-6),sprintf('%02u.txt',iMC)];
            fReplaceText([mcParaPath,filesep,sElastixParamFile{iParams}],sprintf('TransformParameters_Gate%02u_f_prior.txt" "TransformParameters_Gate%02u_f_init.txt',iMC,iMC),'SubTransforms');
            
            % backward
            copyfile([mcParaPath,filesep,sElastixParamFileInit{iParams}(1:end-4),'_inverse.txt'], [mcParaPath,filesep,sElastixParamFileInit{iParams}(1:end-6),sprintf('%02u_inverse.txt',iMC)],'f');
            fReplaceText([mcParaPath,filesep,sElastixParamFile{iParams}(1:end-4),'_inverse.txt'],sprintf('TransformParameters_Gate%02u_b_prior.txt" "TransformParameters_Gate%02u_b_init.txt',iMC,iMC),'SubTransforms');
        end
    end
    
    % forward deformation field
    if(length(sElastixParamFile) == 3) % composed transformation
        if(bInit)
            if(bSurf)
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,sprintf('Gate%02u.mhd', iMC), ' -out ',mcPath,' -p ', mcParaPath,sElastixParamFile{1}, ' -p ', mcParaPath,sElastixParamFile{2},' -p ', mcParaPath,sElastixParamFile{3},' -fp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_fixed.txt',iMC),' -mp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_moved.txt',iMC)]); % with SURF
            else
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,sprintf('Gate%02u.mhd', iMC), ' -out ',mcPath,' -p ', mcParaPath,sElastixParamFile{1}, ' -p ', mcParaPath,sElastixParamFile{2},' -p ', mcParaPath,sElastixParamFile{3}]); % without SURF
            end
        else
            if(bSurf)
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,sprintf('Gate%02u.mhd', iMC), ' -out ',mcPath,' -t0 ',mcPath,filesep, sprintf('TransformParameters_Gate%02u_f.txt', iMC),' -p ', mcParaPath,sElastixParamFile{1}, ' -p ', mcParaPath,filesep,sElastixParamFile{2},' -p ', mcParaPath,filesep,sElastixParamFile{3},' -fp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_fixed.txt',iMC),' -mp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_moved.txt',iMC)]); % with SURF
            else
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,sprintf('Gate%02u.mhd', iMC), ' -out ',mcPath,' -t0 ',mcPath,filesep, sprintf('TransformParameters_Gate%02u_f.txt', iMC),' -p ', mcParaPath,filesep,sElastixParamFile{1}, ' -p ', mcParaPath,filesep,sElastixParamFile{2},' -p ', mcParaPath,filesep,sElastixParamFile{3}]); % without SURF
            end
        end
        copyfile([mcPath,filesep, 'TransformParameters.0.txt'], [mcPath,filesep, sprintf('TransformParameters_Gate%02u_fr.txt', iMC)],'f'); % rigid
        copyfile([mcPath,filesep, 'TransformParameters.1.txt'], [mcPath,filesep, sprintf('TransformParameters_Gate%02u_fa.txt', iMC)],'f'); % affine
        fReplaceText([mcPath,filesep, sprintf('TransformParameters_Gate%02u_fa.txt', iMC)],[mcPath,filesep, sprintf('TransformParameters_Gate%02u_fr.txt', iMC)]);
        copyfile([mcPath,filesep, 'TransformParameters.2.txt'], [mcPath,filesep, sprintf('TransformParameters_Gate%02u_f.txt', iMC)],'f'); % bspline
        fReplaceText([mcPath,filesep, sprintf('TransformParameters_Gate%02u_f.txt', iMC)],[mcPath,filesep, sprintf('TransformParameters_Gate%02u_fa.txt', iMC)]);
    elseif(length(sElastixParamFile) == 1) % fast mode
        if(bInit)
            if(bSurf)
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,sprintf('Gate%02u.mhd', iMC), ' -out ',mcPath,' -p ', mcParaPath,sElastixParamFile{1},' -fp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_fixed.txt',iMC),' -mp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_moved.txt',iMC)]); % with SURF
            else
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,sprintf('Gate%02u.mhd', iMC), ' -out ',mcPath,' -p ', mcParaPath,filesep,sElastixParamFile{1}]); % without SURF
            end
        else
            if(bSurf)
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,sprintf('Gate%02u.mhd', iMC), ' -out ',mcPath,' -t0 ',mcPath,filesep, sprintf('TransformParameters_Gate%02u_f.txt', iMC),' -p ', mcParaPath,filesep,sElastixParamFile{1},' -fp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_fixed.txt',iMC),' -mp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_moved.txt',iMC)]); % with SURF
            else
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,sprintf('Gate%02u.mhd', iMC), ' -out ',mcPath,' -t0 ',mcPath,filesep, sprintf('TransformParameters_Gate%02u_f.txt', iMC),' -p ', mcParaPath,filesep,sElastixParamFile{1}]); % without SURF
            end
        end         
        copyfile([mcPath,filesep, 'TransformParameters.0.txt'], [mcPath,filesep, sprintf('TransformParameters_Gate%02u_f.txt', iMC)],'f'); % bspline
        fReplaceText([mcPath,filesep, sprintf('TransformParameters_Gate%02u_f.txt', iMC)],'NoInitialTransform');
    else
        % NOP
    end
    if(bCalcDF)
        % retrieve deformation field
        [~,~] = system([sElastixPath, filesep, 'transformix.exe -def all -out ',mcPath,' -tp ', mcPath,filesep,sprintf('TransformParameters_Gate%02u_f.txt', iMC)]);

        tmpPath = pwd;
        cd(mcPath);
        SData = read_mhd('deformationField.mhd');
        SDeform(iMC).dFy = SData.datax./SImg.spacing(2); % in [px]
        SDeform(iMC).dFx = SData.datay./SImg.spacing(1); % in [px]
        SDeform(iMC).dFz = SData.dataz./SImg.spacing(3); % in [px]
        cd(tmpPath);
    end
    if(bSurf)
        dispProgress('MC',((iMC-2)*(size(dImg,3)+2)+size(dImg,3)+1)/iTotal);
    else
        dispProgress('MC',(iMC-2)/iTotal);
    end
    
    % backward deformation field
    if(length(sElastixParamFile) == 3)
        if(bCompose) % needed for compose! -> no t0 for this case
            if(bSurf)
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,'Gate01.mhd -out ',mcPath,' -p ', mcParaPath,filesep,sElastixParamFile{3}(1:end-4),'_inverse.txt -fp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_fixed.txt',iMC),' -mp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_moved.txt',iMC)]); % with SURF % TODO: evtl. bei Parafile die SamplingPoint Metric rausnehmen?!
            else
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,'Gate01.mhd -out ',mcPath,' -p ', mcParaPath,filesep,sElastixParamFile{3}(1:end-4),'_inverse.txt']); % without SURF % TODO: evtl. bei Parafile die SamplingPoint Metric rausnehmen?!
            end
        else
            if(bSurf)
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,'Gate01.mhd -out ',mcPath,' -p ', mcParaPath,filesep,sElastixParamFile{3}(1:end-4),'_inverse.txt -t0 ',mcPath,filesep,sprintf('TransformParameters_Gate%02u_f.txt', iMC),' -fp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_fixed.txt',iMC),' -mp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_moved.txt',iMC)]); % with SURF % TODO: evtl. bei Parafile die SamplingPoint Metric rausnehmen?!
            else
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,'Gate01.mhd -out ',mcPath,' -p ', mcParaPath,filesep,sElastixParamFile{3}(1:end-4),'_inverse.txt -t0 ',mcPath,filesep,sprintf('TransformParameters_Gate%02u_f.txt', iMC)]); % without SURF % TODO: evtl. bei Parafile die SamplingPoint Metric rausnehmen?!
            end
        end
        copyfile([mcPath,filesep, 'TransformParameters.0.txt'], [mcPath,filesep, sprintf('TransformParameters_Gate%02u_b.txt', iMC)],'f');
        fReplaceText([mcPath,filesep, sprintf('TransformParameters_Gate%02u_b.txt', iMC)],'NoInitialTransform');
    elseif(length(sElastixParamFile) == 1)
        if(bCompose)
            if(bSurf)
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,'Gate01.mhd -out ',mcPath,' -p ', mcParaPath,filesep,sElastixParamFile{1}(1:end-4),'_inverse.txt -fp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_fixed.txt',iMC),' -mp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_moved.txt',iMC)]); % with SURF % TODO: evtl. bei Parafile die SamplingPoint Metric rausnehmen?!
            else
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,'Gate01.mhd -out ',mcPath,' -p ', mcParaPath,filesep,sElastixParamFile{1}(1:end-4),'_inverse.txt']); % without SURF % TODO: evtl. bei Parafile die SamplingPoint Metric rausnehmen?!
            end
        else
            if(bSurf)
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,'Gate01.mhd -out ',mcPath,' -p ', mcParaPath,filesep,sElastixParamFile{1}(1:end-4),'_inverse.txt -t0 ',mcPath,filesep,sprintf('TransformParameters_Gate%02u_f.txt', iMC),' -fp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_fixed.txt',iMC),' -mp ',mcPath,filesep,sprintf('PointSetGate%02uTo01_moved.txt',iMC)]); % with SURF % TODO: evtl. bei Parafile die SamplingPoint Metric rausnehmen?!
            else
                [~,~] = system([sElastixPath, filesep, 'elastix.exe -f ',mcPath,filesep,'Gate01.mhd -m ',mcPath,filesep,'Gate01.mhd -out ',mcPath,' -p ', mcParaPath,filesep,sElastixParamFile{1}(1:end-4),'_inverse.txt -t0 ',mcPath,filesep,sprintf('TransformParameters_Gate%02u_f.txt', iMC)]); % without SURF % TODO: evtl. bei Parafile die SamplingPoint Metric rausnehmen?!
            end
        end
        copyfile([mcPath,filesep, 'TransformParameters.0.txt'], [mcPath,filesep, sprintf('TransformParameters_Gate%02u_b.txt', iMC)],'f');
        fReplaceText([mcPath,filesep, sprintf('TransformParameters_Gate%02u_b.txt', iMC)],'NoInitialTransform');
    else
        % NOP
    end
    if(bCalcDF)
        % retrieve deformation field
        [~,~] = system([sElastixPath, filesep, 'transformix.exe -def all -out ',mcPath,' -tp ', mcPath,filesep,sprintf('TransformParameters_Gate%02u_b.txt', iMC)]);
        
        tmpPath = pwd;
        cd(mcPath);
        SData = read_mhd('deformationField.mhd');
        SDeform(iMC).dBy = SData.datax./SImg.spacing(2); % in [px]
        SDeform(iMC).dBx = SData.datay./SImg.spacing(1); % in [px]
        SDeform(iMC).dBz = SData.dataz./SImg.spacing(3); % in [px]
        cd(tmpPath);
    end
    if(bSurf)
        dispProgress('MC',((iMC-2)*(size(dImg,3)+2)+size(dImg,3)+2)/iTotal);        
    else
        dispProgress('MC',((iMC-2)*2)/iTotal);
    end
end
dispProgress('MC','Close');

if(bCalcDF)
   SImg = SDeform; 
end


end