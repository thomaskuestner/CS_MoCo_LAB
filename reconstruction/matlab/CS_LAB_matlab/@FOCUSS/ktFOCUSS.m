function dImg = ktFOCUSS(obj, kSpace)
% 2D FOCUSS algorithm
%
% input:
% obj       CS reconstruction object (holding all parameters)
% kSpace    subsampled kSpace
%
% output:
% dImg      reconstructed channel individual image
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

% 2D, nTime == 1
% here: kSpace: k_y - k_x - cha
% calibSize: k_y - k_x (cell: t)
% fft: f -> t   x -> k_x   y -> k_y

fprintf(' |- prepare k-space and image kernel...\n');

[nPha, nFreq, nCha] = size(kSpace);
if(any(obj.trafo.fftBA(1,:)))
    for iBA=find(obj.trafo.fftBA(1,:))
        kSpace = ifftnshift(kSpace,iBA);
    end
end

% same subsampling for all channels
% lMask = abs(kSpace(:,:,:,1)) > 0; % sparse matrix (CS pattern)
% retrieve mask
lMask = obj.fullMask; % k_y - k_x - cha

% calculate transformation basis
[fTrafo, bTrafo, kernelFTrafo, kernelBTrafo, paraTrafo] = compBasis(kSpace,obj);

% prepare initial estimate
fprintf(' |- prepare initial estimate...\n');
W = kSpace; % k_y-k_x-cha
helper = true(size(W));
m = [nPha, nFreq, nCha];    
% find the k-space center
if(iscell(obj.calibSize))
    s = obj.calibSize{1}; % automatically determined calibration size
else
    s = obj.calibSize; % predefined calibration size
end
idx = cell(1,length(s));
for n=1:length(s)
    if(mod(s(n),2) == 0)
        idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
    else
        idx{n} = floor(m(n)/2)+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2)-1;
    end

    tmp = [idx{n}(1) <= 0, idx{n}(end) > m(n)];
    if(any(tmp))
        if(all(tmp)), error('crop(): Both index out of bounds'); end;
        hShift = [abs(idx{n}(1)) + 1, idx{n}(end) - m(n)];
        op = {idx{n}(1), idx{n}(end); '+', '-'; '<', '>'; m(n) + 1, 0};
        eval(sprintf('if(op{1,~tmp} %s hShift(tmp) %s %d), idx{n} = idx{n} %s hShift(tmp); else idx{n} = idx{n}(idx{n} %s %d);end;',op{2,tmp}, op{3,tmp}, op{4,tmp}, op{2,tmp}, op{3,~tmp}, op{4,~tmp}));
   end
end
helper(idx{1},idx{2},:) = false;

% windowing for cut-out region
if(obj.window.windowingOn)
    window = windowND(obj.window.type,[length(idx{1}) length(idx{2})],obj.window.windowOpt{1},obj.window.windowOpt{2});
    window = repmat(window,[1 1 nCha]);
end
 
W(helper) = 0; % just take the k-space center (y values); all other values are zero
if(exist('window','var'))
    W(~helper) = W(~helper) .* window(:); 
    clear 'window';
end

W = abs(fTrafo(W)).^(obj.p); % rho values -> y-x-cha new_base
if(strcmp(obj.sScaling,'self'))
    for c=1:nCha
        helper = W(:,:,c);
        W(:,:,c) = W(:,:,c)./max(helper(:));
    end
elseif(strcmp(obj.sScaling,'energy'))
    pixEnergy = abs(W).^2;
    totEnergy = squeeze(sum(sum(pixEnergy,1),2));
    clear 'pixEnergy';
    for c=1:nCha
        W(:,:,c) = W(:,:,c)./totEnergy(c);
    end
end
clear 'helper';
if(obj.lambdaCalib > 0)
    % prepare image kernel
    kernelImg = cell2mat(shiftdim(obj.kernelImg,-3)); % y - x - cha - cha - (t)
    padsize = size(kernelImg);
    kernelImgL = complex(zeros(size(kernelImg),obj.measPara.precision),zeros(size(kernelImg),obj.measPara.precision));
    if(isa(W,'TRAFO'))
        kernelImgL = TRAFO(kernelImgL,getMeta(W),paraTrafo,'kernel2D');
    end
    for iCha=1:nCha
        if(obj.trafo.kspaceTrafo)
            kernelImgL(:,:,:,iCha) = kernelFTrafo(fftnshift(kernelImg(:,:,:,iCha),unique([obj.trafo.kspaceTrafoFFTdims, obj.trafo.fftdim]),unique([obj.trafo.kspaceTrafoScrambledims, obj.trafo.scrambledim]))); % y-x-cha-cha new_base
        else
            kernelImgL(:,:,:,iCha) = kernelFTrafo(kernelImg(:,:,:,iCha));
        end
    end
    kernelImg = kernelImgL;
    obj.kernelImg = []; % clear in object
    clear 'kernelImgL';
end

fprintf(' |- start reconstruction...\n');

dispProgress('FOCUSS',0,obj.iNOUTER);
for iI = 1:obj.iNOUTER % outer loop for kt-FOCUSS
    % initial estimates for CG
    q = complex(zeros(size(W),obj.measPara.precision),zeros(size(W),obj.measPara.precision));
    rho = complex(zeros(size(W),obj.measPara.precision),zeros(size(W),obj.measPara.precision));
    d = complex(zeros(size(W),obj.measPara.precision),zeros(size(W),obj.measPara.precision)); % with this intialization => d_0 = - g_0
    g_old = complex(ones (size(W),obj.measPara.precision),ones (size(W),obj.measPara.precision)); % this initialization prevents division by zero
    dispProgress('CG',0,obj.iNINNER);
    for iJ = 1:obj.iNINNER % inner loop for conjugate gradient
        e = kSpace - lMask.*bTrafo(rho); % ^= v - F*rho; kSpace (t-k_y-x-cha); rho (f-y-x-cha); e (t-k_y-x-cha)
        finished = false(1,nCha);
        for c=1:nCha
            helper = e(:,:,c);
            finished(c) = norm(helper(:)) <= obj.epsilon;
        end
        if(all(finished))
            break;
        end
%         G = -conj(W).*fTrafo(e) .* nPha .* nFreq + obj.lambda .*q + obj.lambdaCalib .* kernelMult(obj.lambdaCalib); % weighting and sum (equation 43 in Jung et al paper);
        G = -conj(W).*fTrafo(e) + obj.lambda .*q + obj.lambdaCalib .* kernelMult(obj.lambdaCalib);
        beta = zeros(size(W),obj.measPara.precision);
        for c=1:nCha
            G_helper = G(:,:,c);
            g_old_helper = g_old(:,:,c);
            beta(:,:,c) = (G_helper(:)'*G_helper(:))/(g_old_helper(:)'*g_old_helper(:)); 
        end
        clear 'G_helper' 'g_old_helper';
        d = beta.*d - G;
        g_old = G;
        z = bTrafo(W.*d).*lMask; 
        alpha = zeros(size(W),obj.measPara.precision);
        for c=1:nCha
            z_helper = z(:,:,c);
            e_helper = e(:,:,c);
            alpha(:,:,c) =  abs((z_helper(:)'*e_helper(:))/(z_helper(:)'*z_helper(:))); % abs: to ensure a step into positive direction
        end
        q = q + alpha.*d;
        rho = W.*q;
        dispProgress('CG',iJ/obj.iNINNER);
        dispProgress('FOCUSS',((iI - 1)*obj.iNINNER + iJ)/(obj.iNOUTER*obj.iNINNER));
    end
    dispProgress('CG','Close');
    if(obj.espresso.state)
%     tmp = permute(abs(fft2shift(rho)) > 0,[1 3 2 4]);
        cs.smplPtrn = obj.fullMask(:,:,1); % y-x
        cs.pfn = obj.espresso.pfn;
        rho = bTrafo(rho);
        if(any(obj.trafo.fftBA(1,:))) % go back to kSpace
            for iBA=find(obj.trafo.fftBA(1,:))
                rho = fftnshift(rho,iBA);
            end
        end
        [rho,~,~] = espresso_recon(permute(rho,[3 1 2]), obj.espresso.iter, true, cs); % in: cha-k_y-k_x-k_z, out: y-z-x-cha (cart)
        rho = permute(rho,[2 3 1]);
%         lMask = permute(lMask,[2 3 1]);
        if(obj.trafo.kspaceTrafo)
            rho = fftnshift(rho,unique([obj.trafo.kspaceTrafoFFTdims, obj.trafo.fftdim]),unique([obj.trafo.kspaceTrafoScrambledims, obj.trafo.scrambledim]));
        end
        rho = kernelFTrafo(rho); % y-z-x-cha (new_base)
    end
    W = abs(rho).^(obj.p);
    if(strcmp(obj.sScaling,'self'))
        for c=1:nCha
            helper = W(:,:,c);
            W(:,:,c) = W(:,:,c)./max(helper(:));
        end
    elseif(strcmp(obj.sScaling,'energy'))
        for c=1:nCha
            W(:,:,c) = W(:,:,c)./totEnergy(c);
        end
    end
    clear 'helper';
end
dispProgress('FOCUSS', 'Close');
rho = W.*q;
dImg = kernelBTrafo(rho); % y-x-cha (cart)
if(obj.trafo.kspaceTrafo)
    dImg = ifftnshift(dImg,obj.trafo.kspaceTrafoFFTdims,obj.trafo.kspaceTrafoScrambledims); % 3D: y and z
end
if(any(obj.trafo.fftBA(2,:)))
    for iBA=find(obj.trafo.fftBA(2,:))
        dImg = ifftnshift(dImg,iBA);
    end
end
dImg = shiftdim(mat2cell(dImg, nPha, nFreq, ones(1,nCha)),1);

    function res = kernelMult(lambdaCalib)
        if(lambdaCalib == 0)
            res = 0;
            return;
        end
        res = complex(zeros(size(kernelImg(:,:,:,1)),obj.measPara.precision),zeros(size(kernelImg(:,:,:,1)),obj.measPara.precision));
        % zero-padding just allowed in cartesian image domain
        if(any(size(W) ~= size(res))) % same holds for q (q and W always have the same size)
            padsizeCurr = padsize(1:end-1);
            resize = true;
            Wtmp = kernelBTrafo(W);
            cropsize = size(Wtmp);
            if(obj.trafo.zeroPad)
                Wtmp = zpad(Wtmp,padsizeCurr);
            else
                Wout = complex(zeros(padsizeCurr,obj.measPara.precision),zeros(padsizeCurr,obj.measPara.precision));
                for iCha=1:nCha
                    Wout(:,:,iCha) = crop(imresize(Wtmp(:,:,iCha),padsizeCurr(1:2),para.rescaleInterp),padsizeCurr(1:3));
                end
                Wtmp = Wout;
                clear 'Wout';
            end
            Wtmp = kernelFTrafo(Wtmp);

            qtmp = kernelBTrafo(q);
            if(obj.trafo.zeroPad)
                qtmp = zpad(qtmp,padsize);
            else
                qout = complex(zeros(padsizeCurr,obj.measPara.precision),zeros(padsizeCurr,obj.measPara.precision));
                for iCha=1:nCha
                    qout(:,:,iCha) = crop(imresize(qtmp(:,:,iCha),padsizeCurr(1:2),para.rescaleInterp),padsizeCurr(1:3));
                end
                qtmp = qout;
                clear 'qout';
            end
            qtmp = kernelFTrafo(qtmp);
        else
            resize = false;
            Wtmp = W;
            qtmp = q;
        end
        for c=1:nCha
            kernel = kernelImg(:,:,:,c);
            kMI = (kernel - ones(size(kernel))).*Wtmp; % kernel - 1   TODO: check abs(.) nötig
    %         res(:,:,:,c) = sum(kMI.*kMI .* q,4);
            res(:,:,c) = sum(conj(kMI) .* kMI .* qtmp,3);
        end
        clear 'Wtmp' 'qtmp' 'padsizeCurr' 'kMI' 'kernel';

    %     res = crop(res,size(q));
        if(resize)
            res = kernelBTrafo(res);
            if(obj.trafo.zeroPad)
                res = crop(res,cropsize); 
            else
                resout = complex(zeros(cropsize,obj.measPara.precision),zeros(cropsize,obj.measPara.precision));
                for iCha=1:nCha
                    resout(:,:,iCha) = crop(imresize(res(:,:,iCha), cropsize(1:2), para.rescaleInterp),cropsize(1:3));
                end
                res = resout;
                clear 'resout';
            end
            res = kernelFTrafo(res);
        end
    end
end