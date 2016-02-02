function [outdata, meta] = initTRAFO( matrix, paraTrafo, inmeta, initialize )
%INITTRAFO initialize special transformation structure without computing
%transformation
% matrix        zeros or ones
% paraTrafo     transformation parameters
% inmeta        input meta information
% initialize    'kernel' (cha-cha structure) | 'default' (1-cha structure)
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

cutOut = '';
for i=1:ndims(matrix)
    if(i <= paraTrafo.shape)
        cutOut = [cutOut,':,'];
    else
        cutOut = [cutOut,'1,'];
    end
end
cutOut = cutOut(1:end-1);
inmatrixSize = size(matrix);
if(isa(matrix,'double'))
    precision = 'double';
else
    precision = 'single';
end

if(~isempty(paraTrafo.permRule))
    if(strcmp(initialize,'kernel'))
        tmpPermRule = max(paraTrafo.permRule)+1:ndims(matrix);
        matrix = permute(matrix,[paraTrafo.permRule, tmpPermRule]);
    else
        matrix = permute(matrix,paraTrafo.permRule);
    end
end

switch paraTrafo.trafoType
    case 'surfacelet'
        % TODO: implement pure 2D case!!!!
%         minPadSize = max([max(paraTrafo.size), max(size(matrix))]); % TODO: anpassen für pyr_mode ~= 1, diese Zeile gilt so nur für pyr_mode = 1
        sizes = [paraTrafo.imgsize(1:paraTrafo.shape), size(matrix)];
        
        if(isempty(paraTrafo.padsize))
            if(strcmp(initialize,'kernel'))
                n = ceil(sizes(1:2*paraTrafo.shape)./[paraTrafo.downSamp, paraTrafo.downSamp]);
                minPadSize = max(n.*[paraTrafo.downSamp, paraTrafo.downSamp]);

                outdata = cell(paraTrafo.size(end),paraTrafo.size(end));
                if(matrix(1) == 0)
                    eval(sprintf('matrix = zpad(matrix(%s),minPadSize*ones(1,paraTrafo.shape));', cutOut));
                else
                    eval(sprintf('matrix = opad(matrix(%s),minPadSize*ones(1,paraTrafo.shape));', cutOut)); 
                end
            else
                n = ceil(sizes(1:2*paraTrafo.shape)./[paraTrafo.downSamp, paraTrafo.downSamp]);
                minPadSize = max(n.*[paraTrafo.downSamp, paraTrafo.downSamp]);

                outdata = cell(1,paraTrafo.size(end));
                if(matrix(1) == 0)
                    eval(sprintf('matrix = zpad(matrix(%s),minPadSize*ones(1,paraTrafo.shape));', cutOut));
                else
                    eval(sprintf('matrix = opad(matrix(%s),minPadSize*ones(1,paraTrafo.shape));', cutOut)); 
                end
            end
        else
            if(strcmp(initialize,'kernel'))
                outdata = cell(paraTrafo.size(4),paraTrafo.size(4));
                if(matrix(1) == 0)
                    eval(sprintf('matrix = zpad(matrix(%s),paraTrafo.padsize*ones(1,paraTrafo.shape));', cutOut));
                else
                    eval(sprintf('matrix = opad(matrix(%s),paraTrafo.padsize*ones(1,paraTrafo.shape));', cutOut)); 
                end
            else
                outdata = cell(1,paraTrafo.size(4));
                if(matrix(1) == 0)
                    eval(sprintf('matrix = zpad(matrix(%s),paraTrafo.padsize*ones(1,paraTrafo.shape));', cutOut));
                else
                    eval(sprintf('matrix = opad(matrix(%s),paraTrafo.padsize*ones(1,paraTrafo.shape));', cutOut)); 
                end
            end
        end
        downFactor = 2;
        nDim = ndims(matrix);
        if(nDim > 3)
            error('initTrafo(): Maximal 3D decomposition allowed at the moment!');
        end
        
        if(paraTrafo.shape == 2)
            [tmpMat, ~] = Surfdec(matrix(:,:,1,1,1), paraTrafo.Pyr_mode, paraTrafo.Lev_array, paraTrafo.HGfname, 'bo', paraTrafo.bo, 'msize', paraTrafo.msize, 'beta', paraTrafo.beta, 'lambda', paraTrafo.lambda);
            % flatten all cell matrices
            tmpMat = flattenCellMatrix(tmpMat);
            
            if(strcmp(paraTrafo.dimensionality, '4D'))
                outdata{1} = cell(1,inmatrixSize(4));
                outMat = cell([size(tmpMat),inmatrixSize(3)]);
                for j=1:inmatrixSize(3)
                    outMat(:,:,j) = tmpMat;
                end
                [outdata{1}{:}] = deal(outMat);
            elseif(strcmp(paraTrafo.dimensionality, '3D') || strcmp(paraTrafo.dimensionality, '2Dt'))
                outdata{1} = cell([size(tmpMat),inmatrixSize(3)]);
                for j=1:inmatrixSize(3)
                    outdata{1}(:,:,j) = tmpMat;
                end
            elseif(strcmp(paraTrafo.dimensionality, '2D'))
                outdata{1} = tmpMat;
            end
            
        elseif(paraTrafo.shape == 3)
            % maximal 3D decomposition !!!
            outdata{1} = cell(size(paraTrafo.Lev_array,2) + 1,1);
            [outdata{1}{1:end-1}] = deal(cell(3,1));
            for i=1:size(outdata{1},1)
                if(i == 2)
                    ldown = paraTrafo.Pyr_mode;
                    matrix = matrix(round(1:ldown:end),round(1:ldown:end),round(1:ldown:end));                 
                elseif(i > 2)
                    % downsampling along all dimension for i > 2
                    matrix = matrix(1:downFactor:end,1:downFactor:end,1:downFactor:end);
                end

                if(i ~= size(outdata{1},1))
                    for j=1:nDim %3D
                        outdata{1}{i}{j} = cell(2^(sum(paraTrafo.Lev_array{i}(j,~ismember(paraTrafo.Lev_array{i}(j,:),-1)))),1);
                        if(any(paraTrafo.Lev_array{i}(j,:) ~= 0 & paraTrafo.Lev_array{i}(j,:) ~= -1))
                            % local downsampling for subbands along specific dimension
                            downDim = find(paraTrafo.Lev_array{i}(j,:) ~= 0 & paraTrafo.Lev_array{i}(j,:) ~= -1);
                            downFactorDim = 2.^paraTrafo.Lev_array{i}(j,downDim);
                            idx = cell(1,nDim);
                            [idx{:}] = deal(':');
                            for k=1:length(downDim)
                                idx{downDim(k)} = 1:downFactorDim(k):size(matrix,downDim(k));
                            end
                            matrixLocal = matrix(idx{:});
                            [outdata{1}{i}{j}{:}] = deal(matrixLocal);
                        else
                            [outdata{1}{i}{j}{:}] = deal(matrix);
                        end
                    end
                else
                    outdata{1}{i} = matrix;
                end        
            end

            % flatten all cell matrices
            outdata{1} = flattenCellMatrix(outdata{1});
        end
        
        if(strcmp(initialize,'kernel'))
            [outdata{1,2:end}] = deal(outdata{1});
            [outdata{2:end,:}] = deal(outdata{1});
        else
            [outdata{2:end}] = deal(outdata{1});
        end
        meta = inmeta;
        
    case 'wavelet_mat'
        if(strcmp(initialize,'kernel'))
            outdata = cell(paraTrafo.size(end),paraTrafo.size(end));
            meta = outdata;
        else
            outdata = cell(1,paraTrafo.size(end));
            meta = outdata;    
        end
        
        eval(sprintf('matrix = matrix(%s);', cutOut));
        if(paraTrafo.shape == 1)
            [helper, metaLoc] = wavedec(matrix,paraTrafo.waveletStages,paraTrafo.wFilter,'mode',paraTrafo.extMode);
        elseif(paraTrafo.shape == 2)
            [helper, metaLoc] = wavedec2(matrix,paraTrafo.waveletStages,paraTrafo.wFilter);
            repVec = paraTrafo.imgsize(3);
%             if(paraTrafo.imgsize(4) > 1)
%                 repVec = [repVec, paraTrafo.imgsize(4)];
%             end
            metaLoc = repmat(metaLoc, [1 1 repVec]);
        elseif(paraTrafo.shape == 3)
            helper = wavedec3(matrix,paraTrafo.waveletStages,paraTrafo.wFilter,'mode',paraTrafo.extMode);
            metaLoc = rmfield(helper,'dec');
            helper = helper.dec;
        end
        
        if(strcmp(paraTrafo.dimensionality,'4D')) % t-y-z-x-cha
            if(paraTrafo.shape == 3)
                tmpMat = cell(7*paraTrafo.waveletStages+1,paraTrafo.imgsize(4));
            end
            for i=1:paraTrafo.imgsize(4)
                if(paraTrafo.shape == 3)
                    tmpMat(:,i) = helper;
                    continue;
                elseif(paraTrafo.shape == 2)
%                     tmpMat = zeros(paraTrafo.waveletStages+2,2,paraTrafo.imgsize(3),paraTrafo.imgsize(img4));
                    tmpMat = zeros([size(helper),paraTrafo.imgsize(3),paraTrafo.imgsize(img4)],precision);
                end
                for j=1:paraTrafo.imgsize(3)
                    if(paraTrafo.shape == 2)
                        tmpMat(:,:,j,i) = helper;
                        continue;
                    elseif(paraTrafo.shape == 1)
                        tmpMat = zeros(paraTrafo.waveletStages+2,paraTrafo.imgsize(3),paraTrafo.imgsize(3),paraTrafo.imgsize(img4),precision);
                    end
                    for k=1:paraTrafo.imgsize(2)
                        if(paraTrafo.shape == 1)
                            tmpMat(:,k,j,i) = helper;
                        end
                    end
                end
            end
        elseif(strcmp(paraTrafo.dimensionality,'3D') || strcmp(paraTrafo.dimensionality,'2Dt'))
            if(paraTrafo.shape == 3)
                tmpMat = helper;
            elseif(paraTrafo.shape == 2)
%                 tmpMat = zeros(paraTrafo.waveletStages+2,2,paraTrafo.imgsize(3));
                tmpMat = zeros([size(helper),paraTrafo.imgsize(3)],precision);
            end
            for j=1:paraTrafo.imgsize(3)
                if(paraTrafo.shape == 2)
                    tmpMat(:,:,j) = helper;
                    continue;
                elseif(paraTrafo.shape == 1)
                    tmpMat = zeros(paraTrafo.waveletStages+2,paraTrafo.imgsize(2),paraTrafo.imgsize(3),precision);
                end
                for k=1:paraTrafo.imgsize(2)
                    if(paraTrafo.shape == 1)
                        tmpMat(:,k,j) = helper;
                    end
                end
            end
        elseif(strcmp(paraTrafo.dimensionality,'2D'))
            if(paraTrafo.shape == 2)
                tmpMat = helper;
            elseif(paraTrafo.shape == 1)
                tmpMat = zeros(paraTrafo.waveletStages+2,paraTrafo.imgsize(2),precision);
            end
            for j=1:paraTrafo.imgsize(2)
                if(paraTrafo.shape == 1)
                    tmpMat(:,j) = helper;
                end
            end
        end           
        
        [outdata{:}] = deal(tmpMat); 
        [meta{:}] = deal(metaLoc);
%         [meta{:}] = deal(rmfield(helper,'dec'));
        
    case 'curvelab'
        if(strcmp(initialize,'kernel'))
            outdata = cell(paraTrafo.imgsize(end),paraTrafo.imgsize(end));
%             matrix = matrix(:,:,:,1,1);
        else
            outdata = cell(1,paraTrafo.imgsize(end));
%             matrix = matrix(:,:,:,1);
        end
       
        %%% TODO: construct matrix directly -> don't transform it, may lead
        %%% to errors -> we must have a one/zero everywhere in the final
        %%% output structure
        eval(sprintf('matrix = matrix(%s);', cutOut));
        if(paraTrafo.shape == 3)
            helper = fdct3d_forward_mex(size(matrix, 1), size(matrix, 2), size(matrix, 3), paraTrafo.nbscales, paraTrafo.nbdstz_coarse, paraTrafo.allCurvelets, matrix);
        elseif(paraTrafo.shape == 2)
            helper = fdct_usfft_mex(size(matrix,1), size(matrix,2), paraTrafo.nbscales, paraTrafo.nbdstz_coarse, paraTrafo.allCurvelets, matrix);
        end
        
        tmpMat = cell(1, paraTrafo.imgsize(end));
        if(strcmp(paraTrafo.dimensionality,'4D'))
            if(paraTrafo.shape == 3)
                tmpMat = cell(paraTrafo.imgsize(4),1);
            elseif(paraTrafo.shape == 2)
                tmpMat = cell(paraTrafo.imgsize(4),paraTrafo.imgsize(3));
            end
            for i=1:paraTrafo.imgsize(4)
                if(paraTrafo.shape == 3)
                    tmpMat{i} = helper;
                    continue;
                end
                for j=1:paraTrafo.imgsize(3)
                    if(paraTrafo.shape == 2)
                        tmpMat{i,j} = helper;
                    end
                end
            end
        elseif(strcmp(paraTrafo.dimensionality,'3D') || strcmp(paraTrafo.dimensionality,'2Dt'))
            if(paraTrafo.shape == 3)
                tmpMat = helper;
            elseif(paraTrafo.shape == 2)
                tmpMat = cell(paraTrafo.imgsize(3),1);
                for j=1:paraTrafo.imgsize(3)
                    tmpMat{j} = helper;
                end
            end
        elseif(strcmp(paraTrafo.dimensionality,'2D'))
            if(paraTrafo.shape == 2)
                tmpMat = helper;
            end
        end           
              
        [outdata{:}] = deal(tmpMat); 
        meta         = inmeta; 
        
	case 'mellin'
        % happens only for paraTrafo.shape == 2
        if(strcmp(initialize,'kernel'))
            outdata = cell(paraTrafo.imgsize(end),paraTrafo.imgsize(end));
%             matrix = matrix(:,:,1,1,1);
        else
            outdata = cell(1,paraTrafo.imgsize(end));
%             matrix = matrix(:,:,1,1);
        end
       
        %%% TODO: construct matrix directly -> don't transform it, may lead
        %%% to errors -> we must have a one/zero everywhere in the final
        %%% output structure
        eval(sprintf('matrix = matrix(%s);', cutOut));
        if(paraTrafo.shape == 2)
            dim = [size(matrix,2)*paraTrafo.trafoSize(2), size(matrix,1)*paraTrafo.trafoSize(1)];
            helper = AFMT2(matrix, paraTrafo.sigma, dim(1), dim(2), paraTrafo.extrapVal, paraTrafo.interp, 'F-AFMT');
        else
            error('initTrafo(): Unknown Mellin initialization');
        end
       
        if(strcmp(paraTrafo.dimensionality,'4D'))
            helper = repmat(helper,[1 1 paraTrafo.imgsize(3,4)]);
        elseif(strcmp(paraTrafo.dimensionality,'3D') || strcmp(paraTrafo.dimensionality,'2Dt'))
            helper = repmat(helper,[1 1 paraTrafo.imgsize(3)]);
%         elseif(strcmp(paraTrafo.dimensionality,'2D'))
            % helper = helper;
        end 
             
        [outdata{:}] = deal(helper); 
        meta         = inmeta; 
        
end

end

