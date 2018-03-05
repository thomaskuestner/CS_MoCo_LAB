function dImg = ktFOCUSS5D(obj, kSpace)
% 5D FOCUSS recon
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

% here: kSpace: t - k_y - k_z - k_x - g - cha
% calibSize: k_y - k_x - k_z (cell: t)
% fft: f -> t   x -> k_x   y -> k_y

fprintf(' |- prepare k-space and image kernel...\n');

[nTime, nPha, nZ, nFreq, nSmp, nCha] = size(kSpace);

% retrieve mask
obj.fullMask = permute(obj.fullMask, [5 1 2 3 6 4]); % t-k_y-k_z-k_x/x-g-cha

% calculate transformation basis
[fTrafoBase, bTrafoBase, kernelFTrafoBase, kernelBTrafoBase, paraTrafo] = compBasis(kSpace,obj);
sigma = cell(1,nCha);
for iCha=1:nCha
    kSpaceIn = permute(ifftnshift(kSpace(:,:,:,:,:,iCha),paraTrafo.fftdim,paraTrafo.scrambledim),[1 2 3 5 4]); % -> t-y-z-g-x
    sigma{iCha} = zeros(size(kSpaceIn,1)*size(kSpaceIn,2)*size(kSpaceIn,3),size(kSpaceIn,4),size(kSpaceIn,5));
    
    for iFreq = 1:nFreq
        sigma{iCha}(:,:,iFreq) = reshape(kSpaceIn(:,:,:,:,iFreq),size(sigma{iCha},1),size(sigma{iCha},2));
    end
end

pc = zeros(nSmp,nSmp,size(sigma{1},3),nCha);
for iCha = 1:nCha
    for iFreq = 1:nFreq
        pc(:,:,iFreq,iCha) = princomp(sigma{iCha}(:,:,iFreq));
    end
end

clear 'sigma' 'kSpaceIn';

% prepare initial estimate
fprintf(' |- prepare initial estimate...\n');
W = kSpace; % t-k_y-k_z-x-g-cha
helper = true(size(W));
m = [nPha, nZ, nFreq, nSmp, nCha];
if(iscell(obj.calibSize))
    % automatically determined calibration size
    loop = nTime;
    pos = 't,idx{1},idx{2},idx{3},:,:';
else
    % predefined calibration size
    loop = 1;
    pos = ':,idx{1},idx{2},idx{3},:,:';
end
for t=1:loop
    % find the k-space center
    if(iscell(obj.calibSize))
        s = obj.calibSize{t};
    else
        s = obj.calibSize;
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
    eval(sprintf('helper(%s) = false;',pos));
end 
kSpaceCenter = ~helper;

% windowing for cut-out region
if(obj.window.windowingOn)
    window = windowND(obj.window.type,[length(idx{1}) length(idx{2}) length(idx{3})],obj.window.windowOpt{1},obj.window.windowOpt{2});
    window = permute(repmat(window,[1 1 1 nTime nSmp nCha]), [4 1 2 3 5 6]);
end

W(helper) = 0; % just take the k-space center (y values); all other values are zero
clear 'helper' 'idx' 'hShift' 'tmp' 'm' 's';
if(exist('window','var'))
    W(kSpaceCenter) = W(kSpaceCenter) .* window(:); 
    clear 'window';
end

W = abs(fTrafo(W)).^(obj.p); % rho values (-> f-y-z-x-cha space)
if(strcmp(obj.sScaling,'self'))
    for c=1:nCha
        helper = W(:,:,:,:,:,c);
        W(:,:,:,:,:,c) = W(:,:,:,:,:,c)./max(helper(:));
    end
elseif(strcmp(obj.sScaling,'energy'))
    pixEnergy = abs(W).^2;
    totEnergy = squeeze(sum(sum(sum(sum(sum(pixEnergy,1),2),3),4),5));
    clear 'pixEnergy';
    for c=1:nCha
        W(:,:,:,:,:,c) = W(:,:,:,:,:,c)./totEnergy(c);
    end
end
clear 'helper';
if(obj.lambdaCalib > 0)
    % prepare image kernel
    kernelImg = zeros(size(obj.kernelImg));
    padsize = size(kernelImg);
    if(isa(W,'TRAFO'))
        kernelImg = TRAFO(kernelImg, getMeta(W), paraTrafo, 'kernel');
    end
    for iCha=1:nCha
        if(obj.trafo.kspaceTrafo)           
            kernelImg(:,:,:,:,:,iCha) = kernelFTrafo(fftnshift(obj.kernelImg(:,:,:,:,:,iCha),unique([obj.trafo.kspaceTrafoFFTdims, obj.trafo.fftdim]),unique([obj.trafo.kspaceTrafoScrambledims, obj.trafo.scrambledim])));
        else
            kernelImg(:,:,:,:,:,iCha) = kernelFTrafo(obj.kernelImg(:,:,:,:,:,iCha)); % => t-y-z-x-cha-cha (new_base)
        end
    end
    obj.kernelImg = []; % clear in object
end

fprintf(' |- start reconstruction...\n');

dispProgress('FOCUSS',0,obj.iNOUTER);
for iI = 1:obj.iNOUTER % outer loop for kt-FOCUSS
    % initial estimates for CG
    q = complex(zeros(size(W),obj.measPara.precision),zeros(size(W),obj.measPara.precision));
    rho = complex(zeros(size(W),obj.measPara.precision),zeros(size(W),obj.measPara.precision));
    d = complex(zeros(size(W),obj.measPara.precision),zeros(size(W),obj.measPara.precision)); % with this intialization => d_0 = - g_0
    g_old = complex(ones(size(W),obj.measPara.precision),ones(size(W),obj.measPara.precision)); % this initialization prevents division by zero
    dispProgress('CG',0,obj.iNINNER);
    for iJ = 1:obj.iNINNER % inner loop for conjugate gradient
        e = kSpace - obj.fullMask.*bTrafo(rho); % ^= v - F*rho; kSpace (t-k_y-k_z-x-cha); rho (f-y-z-x-cha); e (t-k_y-k_z-x-g-cha)
        finished = false(1,nCha);
        for c=1:nCha
            helper = e(:,:,:,:,:,c);
            finished(c) = norm(helper(:)) <= obj.epsilon;
        end
        clear 'helper';
        if(all(finished))
            break;
        end
        G = -conj(W).*fTrafo(e) + obj.lambda.*q + obj.lambdaCalib .* kernelMult(obj.lambdaCalib) + obj.lambdaTV .* gradTV(obj.lambdaTV) + obj.lambdaESPReSSo .* gradESPReSSO(obj.lambdaESPReSSo); % weighting and sum (equation 43 in Jung et al paper); +lambda .* q can be added; ".*iNLin.*iNPha" from IFFT prefactors
        beta = zeros(size(W),obj.measPara.precision);
        for c=1:nCha
            G_helper = G(:,:,:,:,:,c);
            g_old_helper = g_old(:,:,:,:,:,c);
            beta(:,:,:,:,:,c) = (G_helper(:)'*G_helper(:))/(g_old_helper(:)'*g_old_helper(:)); 
        end
        if(~isempty(paraTrafo.permRule)), beta = permute(beta,paraTrafo.permRule); end
        clear 'G_helper' 'g_old_helper';
        d = beta.*d - G;
%         d_old = d;
        g_old = G;
        z = bTrafo(W.*d).*obj.fullMask;
        alpha = zeros(size(W),obj.measPara.precision);
        for c=1:nCha
            z_helper = z(:,:,:,:,:,c);
            e_helper = e(:,:,:,:,:,c);
            alpha(:,:,:,:,:,c) =  (z_helper(:)'*e_helper(:))/(z_helper(:)'*z_helper(:));
%             alpha(:,:,:,:,c) =  abs((z_helper(:)'*e_helper(:))/(z_helper(:)'*z_helper(:))); % abs: to ensure a step into positive direction
        end
        if(~isempty(paraTrafo.permRule)), alpha = permute(alpha,paraTrafo.permRule); end
        clear 'z_helper' 'e_helper';
        q = q + alpha.*d;
        rho = W.*q;
        dispProgress('CG',iJ/obj.iNINNER);
        dispProgress('FOCUSS',((iI - 1)*obj.iNINNER + iJ)/(obj.iNOUTER*obj.iNINNER));
    end
    dispProgress('CG','Close');
    if(obj.espresso.state && obj.espresso.reconType == 2) % t-y-z-x       %%%% NOT WORKING
        cs.pfn = obj.espresso.pfn;
        rho = bTrafo(rho);
        if(any(obj.trafo.fftBA(1,:))) % go back to kSpace
            for iBA=find(obj.trafo.fftBA(1,:))
                rho = fftnshift(rho,iBA);
            end
        end
        for t=1:nTime
            cs.smplPtrn = permute(obj.fullMask(t,:,:,:,1),[2 4 3 1]); % -> y-x-z
            rho(t,:,:,:,:) = permute(espresso_recon(permute(rho(t,:,:,:,:),[5 2 4 3 1]), obj.espresso.iter, false, cs),[5 2 4 3 1]); % in: cha-k_y-k_x-k_z, out: y-z-x-cha (cart)
        end
        if(obj.trafo.kspaceTrafo)
            rho = fftnshift(rho,unique([obj.trafo.kspaceTrafoFFTdims, obj.trafo.fftdim]),unique([obj.trafo.kspaceTrafoScrambledims, obj.trafo.scrambledim]));
        end
        rho = kernelFTrafo(rho); % t-y-z-x-cha (new_base)
    end
    W = abs(rho).^(obj.p);
    if(strcmp(obj.sScaling,'self'))
        for c=1:nCha
            helper = W(:,:,:,:,:,c);
            W(:,:,:,:,:,c) = W(:,:,:,:,:,c)./max(helper(:));
        end
    elseif(strcmp(obj.sScaling,'energy'))
        for c=1:nCha
            W(:,:,:,:,:,c) = W(:,:,:,:,:,c)./totEnergy(c);
        end
    end
    clear 'helper';
end
dispProgress('FOCUSS', 'Close');
rho = W.*q;
dImg = kernelBTrafo(rho); % => t-y-z-x-g-cha (cart)
if(obj.trafo.kspaceTrafo)
    dImg = ifftnshift(dImg,obj.trafo.kspaceTrafoFFTdims,obj.trafo.kspaceTrafoScrambledims); % 3D: y and z
end
if(any(obj.trafo.fftBA(2,:)))
    for iBA=find(obj.trafo.fftBA(2,:))
        dImg = ifftnshift(dImg,iBA);
    end
end
if(obj.trafo.trafodim(1,4)) % fft operation along time
    dImg = fftnshift(dImg,1,[]);
end
dImg = permute(dImg, [2 4 3 1 5 6]); % => y-x-z-t-g-cha
dImg = shiftdim(mat2cell(dImg, nPha, nFreq, nZ, nTime, nSmp, ones(1,nCha)),4);

    function img = fTrafo(img)
        % (explicit) PCA along 5th dimension g
        img = ifftnshift(img,paraTrafo.fftdim,paraTrafo.scrambledim);
        imgsize = size(img);
        img = permute(img, [1 2 3 5 4 6]); % -> t-y-z-g-x-cha
        img = permute(reshape(img,[imgsize(1)*imgsize(2)*imgsize(3), imgsize(5), imgsize(4), imgsize(6)]),[2 1 3 4]); % -> g - t*y*z - x - cha
        img = permute(mtimesx(pc,'C',img),[2 1 3 4]);
        img = permute(reshape(img,[imgsize(1), imgsize(2), imgsize(3), imgsize(5), imgsize(4), imgsize(6)]),[1 2 3 5 4 6]); % -> t-y-z-x-g-cha
        img = kernelFTrafoBase(img);
    end

    function img = kernelFTrafo(img)
        imgsize = size(img);
        img = permute(img, [1 2 3 5 4 6]); % -> t-y-z-g-x-cha
        img = permute(reshape(img,[imgsize(1)*imgsize(2)*imgsize(3), imgsize(5), imgsize(4), imgsize(6)]),[2 1 3 4]); % -> g - t*y*z - x - cha
        img = permute(mtimesx(pc,'C',img),[2 1 3 4]);
        img = permute(reshape(img,[imgsize(1), imgsize(2), imgsize(3), imgsize(5), imgsize(4), imgsize(6)]),[1 2 3 5 4 6]); % -> t-y-z-x-g-cha
        img = kernelFTrafoBase(img);
    end

    function img = bTrafo(img)
        img = bTrafoBase(img);
        imgsize = size(img);
        img = permute(img, [1 2 3 5 4 6]); % -> t-y-z-g-x-cha
        img = permute(reshape(img,[imgsize(1)*imgsize(2)*imgsize(3), imgsize(5), imgsize(4), imgsize(6)]),[2 1 3 4]); % -> g - t*y*z - x - cha
        img = permute(mtimesx(pc,img),[2 1 3 4]);
        img = permute(reshape(img,[imgsize(1), imgsize(2), imgsize(3), imgsize(5), imgsize(4), imgsize(6)]),[1 2 3 5 4 6]); % -> t-y-z-x-g-cha
    end

    function img = kernelBTrafo(img)
        img = kernelBTrafoBase(img);
        imgsize = size(img);
        img = permute(img, [1 2 3 5 4 6]); % -> t-y-z-g-x-cha
        img = permute(reshape(img,[imgsize(1)*imgsize(2)*imgsize(3), imgsize(5), imgsize(4), imgsize(6)]),[2 1 3 4]); % -> g - t*y*z - x - cha
        img = permute(mtimesx(pc,img),[2 1 3 4]);
        img = permute(reshape(img,[imgsize(1), imgsize(2), imgsize(3), imgsize(5), imgsize(4), imgsize(6)]),[1 2 3 5 4 6]); % -> t-y-z-x-g-cha
    end

    function res = kernelMult(lambdaCalib) %%%% NOT WORKING => TODO
        if(lambdaCalib == 0)
            res = 0;
            return;
        end
        res = complex(zeros(size(kernelImg(:,:,:,:,:,1)),obj.measPara.precision),zeros(size(kernelImg(:,:,:,:,:,1)),obj.measPara.precision));
        if(any(size(W) ~= size(res))) % same holds for q (q and W always have the same size)
            padsizeCurr = padsize(1:end-1);
            resize = true;
            Wtmp = kernelBTrafo(W);
            cropsize = size(Wtmp);
            if(obj.trafo.zeroPad)
                Wtmp = zpad(Wtmp,padsizeCurr);
            else
                Wout = complex(zeros(padsizeCurr,obj.measPara.precision),zeros(padsizeCurr,obj.measPara.precision));
                RI = imref3d(size(squeeze(Wtmp(1,:,:,:,1))),1,1,1);
                scaleFactor = size(squeeze(Wout(1,:,:,:,1)))./size(squeeze(Wtmp(1,:,:,:,1)));
                tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                for t=1:nTime
                    for iCha=1:nCha
                        Wout(t,:,:,:,iCha) = crop(imwarp(Wtmp(t,:,:,:,iCha),RI,tFormResize,para.rescaleInterp),padsizeCurr(2:4));
                    end
                end
                Wtmp = Wout;
                clear 'Wout';
            end
            Wtmp = kernelFTrafo(Wtmp);

            qtmp = kernelBTrafo(q);
            if(obj.trafo.zeroPad)
                qtmp = zpad(qtmp,padsizeCurr);
            else
                qout = complex(zeros(padsizeCurr,obj.measPara.precision),zeros(padsizeCurr,obj.measPara.precision));
                RI = imref3d(size(squeeze(qtmp(1,:,:,:,1))),1,1,1);
                scaleFactor = size(squeeze(qout(1,:,:,:,1)))./size(squeeze(qtmp(1,:,:,:,1)));
                tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                for t=1:nTime
                    for iCha=1:nCha
                        qout(t,:,:,:,iCha) = crop(imwarp(qtmp(t,:,:,:,iCha),RI,tFormResize,para.rescaleInterp),padsizeCurr(2:4));
                    end
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
            kernel = kernelImg(:,:,:,:,:,c);
            kMI = (kernel - ones(size(kernel),obj.measPara.precision)).*zpad(Wtmp,size(kernel)); % kernel - 1   TODO: check abs(.) nötig
    %         res(:,:,:,c) = sum(kMI.*kMI .* q,4);
            res(:,:,:,:,c) = sum(conj(kMI) .* kMI .* zpad(qtmp,size(kernel)),5);
        end
        clear 'Wtmp' 'qtmp' 'padsizeCurr' 'kMI' 'kernel';
        
        if(resize)
            res = kernelBTrafo(res);
            if(obj.trafo.zeroPad)
                res = crop(res,cropsize); 
            else
                resout = complex(zeros(cropsize,obj.measPara.precision),zeros(cropsize,obj.measPara.precision));
                RI = imref3d(size(squeeze(res(1,:,:,:,1))),1,1,1);
                scaleFactor = size(squeeze(resout(1,:,:,:,1)))./size(squeeze(res(1,:,:,:,1)));
                tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                for t=1:nTime
                    for iCha=1:nCha
                        resout(t,:,:,:,iCha) = crop(imwarp(res(t,:,:,:,iCha),RI,tFormResize,para.rescaleInterp),cropsize(2:4));
                    end
                end
                res = resout;
                clear 'resout';
            end
            res = kernelFTrafo(res);
        end
    end

    function tv = gradTV(lambdaTV) %%%% NOT WORKING => TODO
        % calculate gradient of spatial TV from image
        if(lambdaTV == 0)
            tv = 0;
            return;
        end
        rhotmp = kernelBTrafo(W.*q);
        if(nnz(rhotmp) == 0)
            tv = 0; % scalar needed for TRAFO object addition
            return;
        end
        tv = complex(zeros(size(rhotmp),obj.measPara.precision),zeros(size(rhotmp),obj.measPara.precision));
        type = 2;
        l1Smooth = 1e-15;
        p = 1;

        for c=1:nCha
            for t=1:nTime
                x = rhotmp(:,:,:,t,c);     
                if(type == 1)
                    % formula

                    Dx = x([2:end,end],:,:) - x;
                    Dy = x(:,[2:end,end],:) - x;
                    Dz = x(:,:,[2:end,end]) - x;

                    res = cat(4,Dx,Dy,Dz);
                    Gtv = p*res.*(res.*conj(res) + l1Smooth).^(p/2-1);

                    tmp = Gtv(:,:,:,1);
                    adjDx = tmp([1,1:end-1],:,:) - tmp;
                    adjDx(1,:,:) = -tmp(1,:,:);
                    adjDx(end,:,:) = tmp(end-1,:,:);
                    tmp = Gtv(:,:,:,2);
                    adjDy = tmp(:,[1,1:end-1],:) - tmp;
                    adjDy(:,1,:) = -tmp(:,1,:);
                    adjDy(:,end,:) = tmp(:,end-1,:);
                    tmp = Gtv(:,:,:,3);
                    adjDz = tmp(:,:,[1,1:end-1]) - tmp;
                    adjDz(:,:,1) = -tmp(:,:,1);
                    adjDz(:,:,end) = tmp(:,:,end-1);

                    tv(:,:,:,t,c) = adjDx + adjDy + adjDz;


                elseif(type == 2)
                    % Matlab calculated differentiation
                    % corner points
                    i = 1; j = 1; k = 1;
                    tv(i,j,k,t,c) = -(x(i + 1, j, k) - 3*x(i, j, k) + x(i, j, k + 1) + x(i, j + 1, k))./((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2);
                    i = size(tv,1); j = 1; k = 1;
                    tv(i,j,k,t,c) = (x(i, j, k) - x(i - 1, j, k))./((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2) - (2*x(i, j, k + 1) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2).^(1./2));
                    i = 1; j = size(tv,2); k = 1;
                    tv(i,j,k,t,c) = (x(i, j, k) - x(i, j - 1, k))./((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                    i = size(tv,1); j = size(tv,2); k = 1;
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2).^(1./2));
                    i = 1; j = 1; k = size(tv,3);
                    tv(i,j,k,t,c) = (x(i, j, k) - x(i, j, k - 1))./((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                    i = size(tv,1); j = 1; k = size(tv,3);
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2).^(1./2));
                    i = 1; j = size(tv,2); k = size(tv,3);
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i + 1, j, k))./(2*((x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2));
                    i = size(tv,1); j = size(tv,2); k = size(tv,3);
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j, k)).^2).^(1./2));
                    % outer lanes
                    i = 2:size(tv,1)-1; j = 1; k = 1;
                    tv(i,j,k,t,c) = (x(i, j, k) - x(i - 1, j, k))./((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2) - (x(i + 1, j, k) - 3*x(i, j, k) + x(i, j, k + 1) + x(i, j + 1, k))./((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2);
                    i = 1; j = 2:size(tv,2)-1; k = 1;
                    tv(i,j,k,t,c) = (x(i, j, k) - x(i, j - 1, k))./((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2) - (x(i + 1, j, k) - 3*x(i, j, k) + x(i, j, k + 1) + x(i, j + 1, k))./((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2);
                    i = size(tv,1); j = 2:size(tv,2)-1; k = 1;
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (x(i, j, k) - x(i - 1, j, k))./((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2) - (2*x(i, j, k + 1) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2).^(1./2));
                    i = 2:size(tv,1)-1; j = size(tv,2); k = 1;
                    tv(i,j,k,t,c) = (x(i, j, k) - x(i, j - 1, k))./((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                    i = 2:size(tv,1)-1; j = 1; k = size(tv,3);
                    tv(i,j,k,t,c) = (x(i, j, k) - x(i, j, k - 1))./((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                    i = 2:size(tv,1)-1; j = size(tv,2); k = size(tv,3);
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i + 1, j, k))./(2*((x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2));
                    i = 1; j = 2:size(tv,2)-1; k = size(tv,3);
                    tv(i,j,k,t,c) = (x(i, j, k) - x(i, j, k - 1))./((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2));
                    i = size(tv,1); j = 2:size(tv,2)-1; k = size(tv,3);
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2).^(1./2));
                    i = 1; j = 1; k = 2:size(tv,3)-1;
                    tv(i,j,k,t,c) = (x(i, j, k) - x(i, j, k - 1))./((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2) - (x(i + 1, j, k) - 3*x(i, j, k) + x(i, j, k + 1) + x(i, j + 1, k))./((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2);
                    i = 1; j = size(tv,2); k = 2:size(tv,3)-1;
                    tv(i,j,k,t,c) = (x(i, j, k) - x(i, j - 1, k))./((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2));
                    i = size(tv,1); j = 1; k = 2:size(tv,3)-1;
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (x(i, j, k) - x(i - 1, j, k))./((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2) - (2*x(i, j, k + 1) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2).^(1./2));
                    i = size(tv,1); j = size(tv,2); k = 2:size(tv,3)-1;
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2).^(1./2));
                    % outer planes
                    i = 2:size(tv,1)-1; j = 2:size(tv,2)-1; k = 1;
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 6*x(i, j, k) + 2*x(i, j, k + 1) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                    i = 2:size(tv,1)-1; j = 2:size(tv,2)-1; k = size(tv,3);
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2));
                    i = 1; j = 2:size(tv,2)-1; k = 2:size(tv,3)-1;
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 6*x(i, j, k) + 2*x(i, j, k + 1) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                    i = size(tv,1); j = 2:size(tv,2)-1; k = 2:size(tv,3)-1;
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i, j, k + 1) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2));
                    i = 2:size(tv,1)-1; j = 1; k = 2:size(tv,3)-1;
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 6*x(i, j, k) + 2*x(i, j, k + 1) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                    i = 2:size(tv,1)-1; j = size(tv,2); k = 2:size(tv,3)-1;
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2));
                    % inner points
                    i = 2:size(tv,1)-1; j = 2:size(tv,2)-1; k = 2:size(tv,3)-1;
                    tv(i,j,k,t,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 6*x(i, j, k) + 2*x(i, j, k + 1) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));

                end
            end
        end
        tv = kernelFTrafo(tv);
        clear 'Gtv' 'adjDx' 'adjDy' 'adjDz' 'tmp';
    end

    function tv = gradTVTime(lambdaTVTime)
        % calculate the gradient of TV along the time direction
        if(lambdaTVTime == 0)
            tv = 0;
            return;
        end
        rhotmp = kernelBTrafo(W.*q);
        if(nnz(rhotmp) == 0)
            tv = 0; % scalar needed for TRAFO object addition
            return;
        end
        % TODO
    end

    function out = gradESPReSSO(lambdaESPReSSo) %%%% NOT WORKING => TODO
        % calculate ESPReSSo gradient
        if(lambdaESPReSSo == 0 || obj.espresso.reconType ~= 1)
            out = 0;
            return;
        end
        
        rhoTmp = bTrafo(W.*q);
        if(nnz(rhoTmp) == 0)
            rhoTmp = kSpace;
            lInit = true;
        else
            lInit = false;
        end
        if(any(obj.trafo.fftBA(1,:))) % go back to kSpace
            for iBA=find(obj.trafo.fftBA(1,:))
                rhoTmp = fftnshift(rhoTmp,iBA);
            end
        end
        out = complex(zeros(size(rhoTmp),obj.measPara.precision),zeros(size(rhoTmp),obj.measPara.precision));
        for t=1:nTime
            if(~obj.espresso.state)
    %             mask = obj.fullMask(:,:,:,1); % y-z-x
                maskRight = squeeze(obj.fullMask(t,:,:,:,:));
                maskLeft = maskRight;
                maskRight(1:end/2,:,:,:) = false;
                maskLeft(end/2:end,:,:,:) = false;
                maskRight = maskRight | squeeze(kSpaceCenter(t,:,:,:,:)); % y-z-x-cha
                maskLeft = maskLeft | squeeze(kSpaceCenter(t,:,:,:,:));
                cs.smplPtrn = permute(maskRight(:,:,:,1),[1 3 2]);
                cs.pfn = 0.5;
            else      
                cs.smplPtrn = permute(obj.fullMask(t,:,:,:,1),[2 4 3 1]); %tmp(:,:,:,1);
                cs.pfn = obj.espresso.pfn;
            end

            [lmaskSym,pfDim,phase] = espresso_reconGrad(permute(rhoTmp(t,:,:,:,:),[5 2 4 3 1]), cs); % in: cha-k_y-k_x-k_z, out: y-z-x-cha (cart)
            if(lInit)
                phase = zeros(size(phase));
            end

    %         if(obj.trafo.kspaceTrafo)
    %             rhoTmp = fftnshift(rhoTmp,unique([obj.trafo.kspaceTrafoFFTdims, obj.trafo.fftdim]),unique([obj.trafo.kspaceTrafoScrambledims, obj.trafo.scrambledim]));
    %         end
            phase = permute(phase,[1 3 2 4]);
            phase = exp(2i * angle(phase));
            lmaskSym = repmat(squeeze(lmaskSym),[1 1 1 nCha]);
            lmaskSym = permute(lmaskSym,[1 3 2 4]);
            kSpaceL = squeeze(kSpace(t,:,:,:,:));
            if(~obj.espresso.state)
                lmaskConj = maskRight(end:-1:1,:,:,:);

                kSpaceCombi = kSpaceL;
                kSpaceCombi(lmaskConj) = conj(kSpaceL(maskRight));
                lmaskConj = maskLeft(end:-1:1,:,:);
                kSpaceCombi(lmaskConj) = conj(kSpaceL(maskLeft));

            else
                lmaskConj = xor(lmaskSym,squeeze(obj.fullMask(t,:,:,:,:)));

                if(pfDim == 1) %y 
                    lmaskLower = lmaskConj(end:-1:1,:,:,:);
                elseif(pfDim == 2) % x
                    lmaskLower = lmaskConj(:,:,end:-1:1,:);
                elseif(pfDim == 3) % z
                    lmaskLower = lmaskConj(:,end:-1:1,:,:);
                end

                kSpaceCombi = kSpaceL;
                kSpaceCombi(lmaskConj) = conj(kSpaceL(lmaskLower));
            end

%             out(t,:,:,:,:) = 0.5 * (conj(W) .* q .* (1+phase) + conj(W) .* W .* conj(q) .* phase - fTrafo(kSpaceCombi) .* (1+phase));     
            out(t,:,:,:,:) = 2 .* conj(q) + conj(W) .* W .* q .* conj(phase) - conj(W) .* fTrafo(kSpaceCombi) .* conj(phase);
        end
    end
end