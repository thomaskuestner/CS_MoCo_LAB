function image = postproc_main(obj, imageCha, postproc)
%POSTPROC_MAIN main postprocessing task
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

if(strcmp(obj.measPara.dimension,'2D'))
    nSlices = size(imageCha,1);
    if(nSlices == 1)
        % one 2D slice was acquired
        obj.currSlice = 1;
        if(strcmp(postproc.type,'rss'))
            dimTurn = unique(cellfun(@(x) length(size(x)), imageCha));
            image = rssPost(imageCha, dimTurn);
        elseif(strcmp(postproc.type,'SENSE'))
            image = SENSE(obj,imageCha,postproc);
        end
                    
    else
        % several 2D slices were acquired
        image = cell(1,nSlices);
        dispProgress('Slices', 0, nSlices);           
        for iSli=1:nSlices
            if(all(cellfun(@isempty, imageCha(iSli,:))))
                continue;
            end
            obj.currSlice = iSli;
            
            if(strcmp(postproc.type,'rss'))
                dimTurn = unique(cellfun(@(x) length(size(x)), imageCha(iSli,:)));
                image{iSli} = rssPost(imageCha(iSli,:),dimTurn);
            elseif(strcmp(postproc.type,'SENSE'))
                image{iSli} = SENSE(obj,imageCha(iSli,:),postproc);
            end
            dispProgress('Slices', iSli/nSlices);
        end
        dispProgress('Slices', 'Close');
        image = cell2mat(shiftdim(image,-1));
    end
    
elseif(strcmp(obj.measPara.dimension,'3D'))
    if(strcmp(postproc.type,'rss'))
        image = rssPost(imageCha,3);
    elseif(strcmp(postproc.type,'SENSE'))
        image = SENSE(obj, imageCha, postproc);
    end  

elseif(strcmp(obj.measPara.dimension,'4D'))
    if(strcmp(postproc.type,'rss'))
        image = rssPost(imageCha,4);
    elseif(strcmp(postproc.type,'SENSE'))
        image = SENSE(obj, imageCha, postproc);
    end 
    
elseif(strcmp(obj.measPara.dimension,'5D'))
    if(strcmp(postproc.type,'rss'))
        image = rssPost(imageCha,5);
    elseif(strcmp(postproc.type,'SENSE'))
        % not available yet
    end 
%     image = image(:,:,:,:,1);
    
end

if(strcmp(obj.measPara.dimension,'2D') && obj.measPara.dim(4) > 1)
    image = permute(image,[1 2 4 3]);
end

% crop image according to frequency/phase/slice oversampling, flip and
% correct anisotropic pixel aspect ratio
image = cropPost(image, obj.measPara.oversampling, postproc.aniso, postproc.interp, postproc.turnImage, postproc.FreqOversamplingCorr, obj.measPara);



end

