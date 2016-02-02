function imageCha = recon_main(obj, iRep, iAvg)
% main reconstruction function
% prepare constraints and select algorithm depending on image dimensionality
%
% input:
% obj       CS reconstruction object (holding all parameters)
% iRep      current repetition
% iAvg      current average
%
% output:
% imageCha  individual reconstructed channel images
%
% (c) Thomas Kuestner 
% -------------------------------------------------------------------------

nCha = obj.measPara.dim(5);
nSlices = size(obj.kSpace,1);
nTime = obj.measPara.dim(4);
imageCha = cell(nSlices,nCha);

if(strcmp(obj.measPara.dimension,'2D'))
    % 2D case
    % extract sampling mask (same for all channels, but different for each
    % slice)
    obj.fullMask = cellfun(@(x) abs(x) > 0, obj.kSpace, 'UniformOutput', false); % k_y-k_x-t-cha or k_y-k_x-cha
    % correct for antiAliasing
    if(obj.measPara.oversampling{2,1})
        obj.fullMask = cellfun(@(x) x(:,obj.measPara.oversampling{1,1},:), obj.fullMask, 'UniformOutput', false);
        obj.kSpace = cellfun(@(x) ifftnshift(x,2), obj.kSpace, 'UniformOutput', false); % k_y - x -t or k_y - x
        obj.kSpace = cellfun(@(x) x(:,obj.measPara.oversampling{1,1},:), obj.kSpace, 'UniformOutput', false);
        obj.kSpace = cellfun(@(x) fftnshift(x,2), obj.kSpace, 'UniformOutput', false); % k_y-k_x-t or k_y-k_x
        obj.kSpace = cellfun(@(x,y) x.*y, kSpace, obj.fullMask, 'UniformOutput', false); % just take acquired points
        obj.measPara.dim(2) = size(obj.kSpace{1,1},2);
    end
    
    if(~isreal(kSpace))
        % boolean variable for complex images
        bComplex = true;
        N = 2*obj.measPara.dim(1)*obj.measPara.dim(2);
        fct = sqrt(2)/sqrt(N/2);
    else
        bComplex = false;
        N = obj.measPara.dim(1)*obj.measPara.dim(2);
        fct = sqrt(2)/sqrt(N);
    end
    P = 1:N;
    
    
    dispProgress('Slices', 0, nSlices);
    for iSli = 1:nSlices
        dispProgress('Channels', 0, nCha);
        for iCha=1:nCha
            dispProgress('Time', 0, nTime);
            imageCha{iSli,iCha} = zeros(obj.measPara.dim(1),obj.measPara.dim(2),obj.measPara.dim(4));
            for t=1:nTime
                mask = logical(obj.fullMask{1}(:,:,t)); % same for each channel
                A = @(z) A_fft2(z, mask, P, obj.measPara.dim(1), obj.measPara.dim(2), bComplex);
                At = @(z) At_fft2(z, N, mask, P, obj.measPara.dim(1), obj.measPara.dim(2), bComplex);
                
                kSpaceIn = double(obj.kSpace{iSli,iCha,iRep,iAvg}(:,:,t));
                b = fct .* [real(kSpaceIn(mask)); imag(kSpaceIn(mask))];
                x0 = At(b);
                
                if(strcmp(obj.solver,'TV'))    
                    img = tvqc_logbarrier(x0, A, At, b, obj.epsilon, obj.lbtol, obj.mu, obj.cgtol, obj.cgmaxiter, bComplex, [obj.measPara.dim(1),obj.measPara.dim(2)]);
                elseif(strcmp(obj.solver,'L1'))
                    img = l1qc_logbarrier(x0, A, At, b, obj.epsilon, obj.lbtol, obj.mu, obj.cgtol, obj.cgmaxiter);
                elseif(strcmp(obj.solver,'TVDantzig'))
                    img = tvdantzig_logbarrier(x0, A, At, b, obj.epsilon, obj.lbtol, obj.mu, obj.cgtol, obj.cgmaxiter, bComplex, [obj.measPara.dim(1),obj.measPara.dim(2)]);  
                elseif(strcmp(obj.solver,'L1Dantzig'))
                    img = l1dantzig_pd(x0, A, At, b, obj.epsilon, obj.lbtol, obj.pdmaxiter, obj.cgtol, obj.cgmaxiter);
                end
                    
                if(bComplex)
                    imageCha{iSli,iCha}(:,:,t) = reshape(img(1:N/2) + 1i * img(N/2+1:N),obj.measPara.dim(1),obj.measPara.dim(2));
                else
                    imageCha{iSli,iCha}(:,:,t) = reshape(img,obj.measPara.dim(1),obj.measPara.dim(2));
                end
                dispProgress('Time', t/nTime);
            end
            dispProgress('Time', 'Close');            
            dispProgress('Channels', iCha/nCha);
        end
        dispProgress('Channels', 'Close');
        dispProgress('Slices', iSli/nSlices); 
    end
    dispProgress('Slices', 'Close');
        
elseif(strcmp(obj.measPara.dimension,'3D'))
    % 3D case
    
    % extract sampling mask (same for all channels, but different for each
    % slice)
    obj.fullMask = cellfun(@(x) abs(x) > 0, obj.kSpace, 'UniformOutput', false); % k_y-k_x-k_z (cell: cha)
    % correct for antiAliasing
    if(obj.measPara.oversampling{2,1})
        obj.fullMask = cellfun(@(x) x(:,obj.measPara.oversampling{1,1},:), obj.fullMask, 'UniformOutput', false);
        obj.kSpace = cellfun(@(x) ifftnshift(x,2), obj.kSpace, 'UniformOutput', false); % k_y - x - k_z (cell: cha)
        obj.kSpace = cellfun(@(x) x(:,obj.measPara.oversampling{1,1},:), obj.kSpace, 'UniformOutput', false);
        obj.kSpace = cellfun(@(x) fftnshift(x,2), obj.kSpace, 'UniformOutput', false); % k_y-k_x-k_z
        obj.kSpace = cellfun(@(x,y) x.*y, obj.kSpace, obj.fullMask, 'UniformOutput', false); % just take acquired points
        obj.measPara.dim(2) = size(obj.kSpace{1,1},2);
    end
    
    if(~isreal(obj.kSpace{1}))
        % boolean variable for complex images
        bComplex = true;
        N = 2*obj.measPara.dim(1)*obj.measPara.dim(2)*obj.measPara.dim(3);
        fct = sqrt(2)/sqrt(N/2);
    else
        bComplex = false;
        N = obj.measPara.dim(1)*obj.measPara.dim(2)*obj.measPara.dim(3);
        fct = sqrt(2)/sqrt(N);
    end
    P = 1:N;
    
    mask = logical(obj.fullMask{1}); % same for each channel
    A = @(z) A_fft3(z, mask, P, obj.measPara.dim(1), obj.measPara.dim(2), obj.measPara.dim(3), bComplex);
    At = @(z) At_fft3(z, N, mask, P, obj.measPara.dim(1), obj.measPara.dim(2), obj.measPara.dim(3), bComplex);
    
    dispProgress('Channels', 0, nCha);
    for iCha=1:nCha              
        kSpaceIn = double(obj.kSpace{1,iCha,iRep,iAvg});
        b = fct .* [real(kSpaceIn(mask)); imag(kSpaceIn(mask))];
        x0 = At(b);
        
        if(strcmp(obj.solver,'TV'))    
            img = tvqc_logbarrier(x0, A, At, b, obj.epsilon, obj.lbtol, obj.mu, obj.cgtol, obj.cgmaxiter, bComplex, [obj.measPara.dim(1),obj.measPara.dim(2),obj.measPara.dim(3)]);
        elseif(strcmp(obj.solver,'L1'))
            img = l1qc_logbarrier(x0, A, At, b, obj.epsilon, obj.lbtol, obj.mu, obj.cgtol, obj.cgmaxiter);
        elseif(strcmp(obj.solver,'TVDantzig'))
            img = tvdantzig_logbarrier(x0, A, At, b, obj.epsilon, obj.lbtol, obj.mu, obj.cgtol, obj.cgmaxiter, bComplex, [obj.measPara.dim(1),obj.measPara.dim(2),obj.measPara.dim(3)]);         
        elseif(strcmp(obj.solver,'L1Dantzig'))
            img = l1dantzig_pd(x0, A, At, b, obj.epsilon, obj.lbtol, obj.pdmaxiter, obj.cgtol, obj.cgmaxiter);
        end
                
        if(bComplex)
            imageCha{1,iCha} = reshape(img(1:N/2) + 1i * img(N/2+1:N),[obj.measPara.dim(1),obj.measPara.dim(2),obj.measPara.dim(3)]);
        else
            imageCha{1,iCha} = reshape(img,[obj.measPara.dim(1),obj.measPara.dim(2),obj.measPara.dim(3)]);
        end
        dispProgress('Channels', iCha/nCha);
    end
    dispProgress('Channels', 'Close');
    
elseif(strcmp(obj.measPara.dimension,'4D'))
    % 4D case

end

