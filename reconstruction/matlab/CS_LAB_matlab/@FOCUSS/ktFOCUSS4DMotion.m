function dImg = ktFOCUSS4DMotion(obj, kSpace)
% 4D FOCUSS algorithm with motion field estimation

% here: kSpace: k_y - k_z - k_x - cha - t
% calibSize: k_y - k_x - k_z (cell: t)
% fft: f -> t   x -> k_x   y -> k_y

fprintf(' |- prepare k-space and image kernel...\n');

[nTime, nPha, nZ, nFreq, nCha] = size(kSpace);

% retrieve mask
obj.fullMask = permute(obj.fullMask, [5 1 2 3 4]); % t-k_y-k_z-k_x/x-cha

% correct calibSize -> better initial registration
[obj.calibSize{2:end}] = deal(obj.calibSize{1});

% prepare initial estimate
fprintf(' |- prepare initial estimate...\n');
W = kSpace; % t-k_y-k_z-x-cha
helper = true(size(W));
m = [nPha, nZ, nFreq, nCha];
if(iscell(obj.calibSize))
    % automatically determined calibration size
    loop = nTime;
    pos = 't,idx{1},idx{2},idx{3},:';
else
    % predefined calibration size
    loop = 1;
    pos = ':,idx{1},idx{2},idx{3},:';
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
    window = permute(repmat(window,[1 1 1 nTime nCha]), [4 1 2 3 5]);
end

W(helper) = 0; % just take the k-space center (y values); all other values are zero
clear 'helper' 'idx' 'hShift' 'tmp' 'm' 's';
if(exist('window','var'))
    W(kSpaceCenter) = W(kSpaceCenter) .* window(:);
    clear 'window';
end

% initialize motion correction
if(obj.lambdaMC > 0)
    % sElastixParamFile = {'rigid.txt', 'affine.txt', 'bSpline.txt'};
    if(strcmp(obj.mc.sElastixParams,'normal'))
        sElastixParamFile = {'default.txt'};
    % sElastixParamFileFast = {'bSpline_fast.txt'};
        sElastixParamFileFast = {'default_fast.txt'};
        sElastixParamFileCompose = {'default_composeXX.txt'};
    elseif(strcmp(mc.sElastixParams,'concatenated'))
        sElastixParamFile = {'rigid.txt', 'affine.txt', 'bSpline.txt'};
        sElastixParamFileFast = {'bSpline_fast.txt'};
        sElastixParamFileCompose = {'bSpline_composeXX.txt'};
    end
    if(isempty(obj.mc.savePath))
        obj.mc.savePath = cd;
    end
    if(~exist([obj.mc.savePath,filesep,'elastix'], 'dir'))
        mkdir([obj.mc.savePath,filesep,'elastix',filesep,'para']);
    end
    mcPath = [obj.mc.savePath,filesep,'elastix'];
    currpath = fileparts(fileparts(mfilename('fullpath')));
    mcParaPath = [currpath,filesep,'utils',filesep,'utils_elastix'];
    if(obj.mc.lSURF)% with SURF
        copyfile([mcParaPath,filesep,'SURF'],[mcPath,filesep,'para',filesep,'SURF'],'f'); % copy para file to data directory without blank spaces
        mcParaPath = [mcPath,filesep,'para',filesep,'SURF'];
    else % without SURF
        copyfile([mcParaPath,filesep,'normal'],[mcPath,filesep,'para',filesep,'normal'],'f'); % copy para file to data directory without blank spaces
        mcParaPath = [mcPath,filesep,'para',filesep,'normal'];
    end

    dImg = W; % t-k_y-k_z-x-cha
    iTrafodim = find(obj.trafo.fftBA(1,:) == false);
    iTrafodim = iTrafodim(~ismember(iTrafodim,1)); % except for t
    for iBA=iTrafodim
        dImg = ifftnshift(dImg,iBA);
    end
    clear 'iTrafodim' 'iBA';

    dImgCenter = dImg;
    dImg = sqrt(sum(abs(dImg).^2,5));
    dImg = fTurnImg(dImg,0);
%     dImg = sqrt(sum(abs(permute(dImg,[2 4 3 1 5])).^2,5)); % TODO: Or multichannel input -> elastix: multi-image, parameter:
    % calculate sensemap
    % a) sum(Multi-Cha-Img .* conj(dSensemap),5)
    % b) repmat(Single-Cha-Img,[1 1 1 1 nCha]) .* dSensemap
    if(obj.mc.lSenseMC)
%         dImgCenter = permute(dImgCenter, [4, 2, 3, 1, 5]);
%         dImgCenter = flipdim(dImgCenter,1);
        dSensemap = dImgCenter./(repmat(sqrt(sum(abs(dImgCenter).^2,5)),[1 1 1 1 size(dImgCenter,5)])); % t-y-z-x-cha
        clear 'dImgCenter';
        dWindow = windowND(@hamming, [size(dSensemap,2), size(dSensemap,3), size(dSensemap,4)]);
        dSensemap = fftnshift(dSensemap,2:4);
        dSensemap = dSensemap .* repmat(shiftdim(dWindow,-1), [size(dSensemap,1) 1 1 1 size(dSensemap,5)]);
        dSensemap = ifftnshift(dSensemap,2:4);
        dSensemapFctr = sqrt(sum(conj(dSensemap).*dSensemap,5)).^(-1);
        dSensemapFctr(isinf(dSensemapFctr)) = 0;
        dSensemapMC = repmat(dSensemapFctr,[1 1 1 1 size(dSensemap,5)]) .* dSensemap;
        clear 'dSensemapFctr' 'dWindow' 'dSensemap';
    else
        clear 'dImgCenter';
    end

    SImg.size = [size(dImg, 1), size(dImg, 2), size(dImg, 3)];
    SImg.orientation = eye(3);
    SImg.spacing = [obj.measPara.FOVread/obj.measPara.baseResolution, obj.measPara.FOVphase/obj.measPara.lines, obj.measPara.sliceThickness];
    SImg.origin = [0; 0; 0];

    % perform initial MC estimation
    fprintf(' |- prepare MC estimation\n');
    if(exist([mcPath,filesep,'init'],'dir'))
        copyfile([mcPath,filesep,'init',filesep,'*'],mcPath,'f');
    else
        SImg = fMCReconMain(dImg, SImg, obj.mc.sElastixPath, mcPath, mcParaPath, sElastixParamFile, 1);
        % manually copy files
        % fCreateLoopImages(dImg,SImg);
    end
    clear 'dImg';
end

if(obj.lSense)
    dImg = W; % t-k_y-k_z-x-cha
    iTrafodim = find(obj.trafo.fftBA(1,:) == false);
    iTrafodim = iTrafodim(ismember(iTrafodim,[2 3 4])); % except for t
    for iBA=iTrafodim
        dImg = ifftnshift(dImg,iBA);
    end
    clear 'iTrafodim' 'iBA';
    dImgCenter = dImg; % t-y-z-x-cha
    dImg = sqrt(sum(abs(dImg).^2,5)); % t-y-z-x
    dSensemap = dImgCenter./(repmat(dImg,[1 1 1 1 size(dImgCenter,5)]));
    clear 'dImgCenter' 'dImg';
    dWindow = windowND(@hamming, [size(dSensemap,2), size(dSensemap,3), size(dSensemap,4)]);
    dSensemap = fftnshift(dSensemap,1:3);
    dSensemap = dSensemap .* permute(repmat(dWindow, [1 1 1 size(dSensemap,1) size(dSensemap,5)]),[4 1 2 3 5]);
    dSensemap = ifftnshift(dSensemap,2:4);
    dSensemapFctr = sqrt(sum(conj(dSensemap).*dSensemap,5)).^(-1);
    dSensemapFctr(isinf(dSensemapFctr)) = 0;
    dSensemap = conj(repmat(dSensemapFctr,[1 1 1 1 size(dSensemap,5)]) .* dSensemap);
    clear 'dSensemapFctr' 'dWindow';
else
    dSensemap = 1;
end

% calculate transformation basis
[fTrafo, bTrafo, kernelFTrafo, kernelBTrafo, paraTrafo] = compBasis(kSpace,obj,dSensemap);

% save initial phase of W for correct back transformation of image in time
% direction, abs(.) enforces smoothing out/projection operation along time
W = fTrafo(W);
% Wphase = angle(W);
W = abs(W).^(obj.p); % rho values (-> f-y-z-x-cha space)
% W = W.^(obj.p);
if(strcmp(obj.sScaling,'self'))
    for c=1:nCha
        helper = W(:,:,:,:,c);
        W(:,:,:,:,c) = W(:,:,:,:,c)./max(helper(:));
    end
elseif(strcmp(obj.sScaling,'energy'))
    pixEnergy = abs(W).^2;
    totEnergy = squeeze(sum(sum(sum(sum(pixEnergy,1),2),3),4));
    clear 'pixEnergy';
    for c=1:nCha
        W(:,:,:,:,c) = W(:,:,:,:,c)./totEnergy(c);
    end
end
clear 'helper';
if(obj.lambdaCalib > 0)
    % prepare image kernel
    kernelImg = zeros(size(obj.kernelImg),obj.measPara.precision);
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
    rho = complex(zeros(size(W),obj.measPara.precision),zeros(size(W),obj.measPara.precision)); clear 'rhoMC';
    d = complex(zeros(size(W),obj.measPara.precision),zeros(size(W),obj.measPara.precision)); % with this intialization => d_0 = - g_0
    g_old = complex(ones(size(W),obj.measPara.precision),ones(size(W),obj.measPara.precision)); % this initialization prevents division by zero
    dispProgress('CG',0,obj.iNINNER);
    for iJ = 1:obj.iNINNER % inner loop for conjugate gradient
        e = kSpace - obj.fullMask.*bTrafo(rho); % ^= v - F*rho; kSpace (t-k_y-k_z-x-cha); rho (f-y-z-x-cha); e (t-k_y-k_z-x-cha)
        finished = false(1,nCha);
        for c=1:nCha
            helper = e(:,:,:,:,c);
            finished(c) = norm(helper(:)) <= obj.epsilon;
        end
        clear 'helper';
        if(all(finished))
            break;
        end
        G = -conj(W).*fTrafo(e) + obj.lambda.*q + obj.lambdaCalib .* kernelMult(obj.lambdaCalib) + obj.lambdaTV .* gradTV(obj.lambdaTV) + obj.lambdaESPReSSo .* gradESPReSSO(obj.lambdaESPReSSo) + obj.lambdaMC .* gradMC(obj.lambdaMC); % weighting and sum (equation 43 in Jung et al paper);
        beta = zeros(size(W),obj.measPara.precision);
        for c=1:nCha
            G_helper = G(:,:,:,:,c);
            g_old_helper = g_old(:,:,:,:,c);
            beta(:,:,:,:,c) = (G_helper(:)'*G_helper(:))/(g_old_helper(:)'*g_old_helper(:));
        end
        if(~isempty(paraTrafo.permRule)), beta = permute(beta,paraTrafo.permRule); end
        clear 'G_helper' 'g_old_helper';
        d = beta.*d - G;
        g_old = G;
        z = bTrafo(W.*d).*obj.fullMask;
        alpha = zeros(size(W),obj.measPara.precision);
        for c=1:nCha
            z_helper = z(:,:,:,:,c);
            e_helper = e(:,:,:,:,c);
            alpha(:,:,:,:,c) =  (z_helper(:)'*e_helper(:))/(z_helper(:)'*z_helper(:));
%             alpha(:,:,:,:,c) =  abs((z_helper(:)'*e_helper(:))/(z_helper(:)'*z_helper(:))); % abs: to ensure a step into positive direction
        end
        if(~isempty(paraTrafo.permRule)), alpha = permute(alpha,paraTrafo.permRule); end
        clear 'z_helper' 'e_helper';
        q = q + alpha.*d;
        rho = W.*q;
        if(obj.lambdaMC > 0)
            % update motion estimation
            % go back to image (cart) => t-y-z-x-cha (cart)
            rhoMC = kernelBTrafo(rho);
            if(obj.trafo.kspaceTrafo)
                rhoMC = ifftnshift(rhoMC,obj.trafo.kspaceTrafoFFTdims,obj.trafo.kspaceTrafoScrambledims); % 3D: y and z
            end
            if(any(obj.trafo.fftBA(2,:)))
                for iBA=find(obj.trafo.fftBA(2,:))
                    rhoMC = ifftnshift(rhoMC,iBA);
                end
            end
            if(obj.trafo.trafodim(1,4)) % fft operation along time
                rhoMC = fftnshift(rhoMC,1,[]);
            end
            rhoMC = sqrt(sum(abs(rhoMC).^2,5)); % TODO: Sensemap?!
            rhoMC = fTurnImg(rhoMC,0);

            SImg = fMCReconMain(rhoMC, SImg, obj.mc.sElastixPath, mcPath, mcParaPath, sElastixParamFileFast);
            % fCreateLoopImages(rhoMC,SImg);
        end
        dispProgress('CG',iJ/obj.iNINNER);
        dispProgress('FOCUSS',((iI - 1)*obj.iNINNER + iJ)/(obj.iNOUTER*obj.iNINNER));
    end
    dispProgress('CG','Close');
    if(obj.espresso.state && obj.espresso.reconType == 2) % t-y-z-x
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
            helper = W(:,:,:,:,c);
            W(:,:,:,:,c) = W(:,:,:,:,c)./max(helper(:));
        end
    elseif(strcmp(obj.sScaling,'energy'))
        for c=1:nCha
            W(:,:,:,:,c) = W(:,:,:,:,c)./totEnergy(c);
        end
    end
    clear 'helper';
end
dispProgress('FOCUSS', 'Close');
rho = W.*q;
dImg = kernelBTrafo(rho); % => f-y-z-x-cha (cart)
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

% for evaluation: M, N, K
obj.measPara.dEval = [prod(size(kSpace)), nnz(obj.fullMask), nnz(q(:) > mean(abs(q(:))) + std(q(:))), nnz(q(:) > mean(abs(q(:))) - std(q(:))), nnz(q(:) > mean(q(:)) + std(q(:))), nnz(q(:) > mean(q(:)) - std(q(:)))];

dImgMC = sqrt(sum(abs(dImg).^2,5));
dImgMC = fTurnImg(dImgMC, 0);
dImg = permute(dImg, [2 4 3 1 5]); % => y-x-z-t-cha
% if(obj.lambdaMC == 0)
%     bInit = true;
% else
%     bInit = false;
% end
% SDeform = fMCReconMain(sqrt(sum(abs(dImg).^25)), SImg, obj.mc.sElastixPath, mcPath, mcParaPath, sElastixParamFile, log2dec(bInit), true); % TODO: Sensemap?!
if(obj.lambdaMC > 0)
    for iMC=2:size(dImg,4)
        copyfile([mcPath,filesep, sprintf('TransformParameters_Gate%02u_f.txt', iMC)],[mcPath,filesep,sprintf('TransformParameters_Gate%02u_f_prior.txt', iMC)]);
        copyfile([mcPath,filesep, sprintf('TransformParameters_Gate%02u_b.txt', iMC)],[mcPath,filesep,sprintf('TransformParameters_Gate%02u_b_prior.txt', iMC)]);
    end
    SImg = fMCReconMain(dImgMC, SImg, obj.mc.sElastixPath, mcPath, mcParaPath, sElastixParamFile, 1, false); % create an "initial" B-spline transformation with normal/SURF param files
    for iMC=2:size(dImg,4)
        copyfile([mcPath,filesep, sprintf('TransformParameters_Gate%02u_f.txt', iMC)],[mcPath,filesep,sprintf('TransformParameters_Gate%02u_f_init.txt', iMC)]);
        copyfile([mcPath,filesep, sprintf('TransformParameters_Gate%02u_b.txt', iMC)],[mcPath,filesep,sprintf('TransformParameters_Gate%02u_b_init.txt', iMC)]);
    end
    % SDeform = fMCReconMain(sqrt(sum(abs(dImg).^2,5)), SImg, obj.mc.sElastixPath, mcPath, [mcPath,filesep,'para',filesep,'compose'], sElastixParamFileCompose, true, true); % compose transformations
    SDeform = fMCReconMain(dImgMC, SImg, obj.mc.sElastixPath, mcPath, mcParaPath, sElastixParamFileCompose, 3, true); % compose transformations

    save([obj.mc.savePath,filesep,'deformdata.mat'],'SDeform');
end

dImg = shiftdim(mat2cell(dImg, nPha, nFreq, nZ, nTime, ones(1,nCha)),3);

    function res = kernelMult(lambdaCalib)
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
            kMI = (kernel - ones(size(kernel),obj.measPara.precision)).*zpad(W,size(kernel)); % kernel - 1   TODO: check abs(.) n�tig
    %         res(:,:,:,c) = sum(kMI.*kMI .* q,4);
            res(:,:,:,:,c) = sum(conj(kMI) .* kMI .* zpad(q,size(kernel)),5);
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

    function tv = gradTV(lambdaTV)
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

    function out = gradESPReSSO(lambdaESPReSSo)
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
            out = 2 .* conj(q) + conj(W) .* W .* q .* conj(phase) - conj(W) .* fTrafo(kSpaceCombi) .* conj(phase);
        end
        clear 'rhoTmp';
    end

    function out = gradMC(lambdaMC)
        % calculate Motion Correction gradient
        if(lambdaMC == 0)
            out = 0;
            return;
        end

        out = complex(zeros(size(W),obj.measPara.precision),zeros(size(W),obj.measPara.precision));
%         rhoTmp = sqrt(sum(abs(permute(rhoTmp,[2 4 3 1 5])).^2,5)); % TODO: Sensemap?!
        tmpPath = pwd;
        cd(mcPath); % due to horrible programming of read_mhd -> will fail otherwise
        if(~exist('rhoMC','var'))
%             rhoRef = permute(out(1,:,:,:,1),[2 4 3 1]); % happens for the first run
            rhoRef = read_mhd('Gate01.mhd'); % happens for the first run
            rhoRef = rhoRef.data;
        else
%             rhoRef = scaleImg(rhoMC(:,:,:,1),[0 2^16-1]);
            rhoRef = rhoMC(:,:,:,1);
        end

        for iMC=2:size(out,1)
            % T_mu
            [~,~] = system([obj.mc.sElastixPath, filesep, 'transformix.exe -in ', [mcPath,filesep,sprintf('Gate%02u.mhd', iMC)],' -out ',mcPath,' -tp ',mcPath,filesep,sprintf('TransformParameters_Gate%02u_f.txt', iMC)]);
            sData = read_mhd('result.mhd');
            % refGate - movedGate
%             SImg.data = int16(rhoRef - sData.data);
            SImg.data = rhoRef - sData.data;
            write_mhd([mcPath,filesep,sprintf('GateDiff%02u.mhd', iMC)], SImg, 'ElementType', 'single');
            % T_mu^H
            [~,~] = system([obj.mc.sElastixPath, filesep, 'transformix.exe -in ', [mcPath,filesep,sprintf('GateDiff%02u.mhd', iMC)],' -out ',mcPath,' -tp ',mcPath,filesep,sprintf('TransformParameters_Gate%02u_b.txt', iMC)]);
            sData = read_mhd('result.mhd');
%             sData.data = sData.data./max(sData.data(:));
            if(obj.mc.lSenseMC)
                out(iMC,:,:,:,:) = repmat(fTurnImg(sData.data,1),[1 1 1 size(out,5)]) .* squeeze(dSensemapMC(iMC,:,:,:,:)); % t-y-z-x-cha
            else
                out(iMC,:,:,:,:) = repmat(fTurnImg(sData.data,1),[1 1 1 size(out,5)]);
            end
%             out(iMC,:,:,:,:) = scaleImg(out(iMC,:,:,:,:), [min(abs(e(:))), max(abs(e(:)))]);
%             out(iMC,:,:,:,:) = scaleImg(out(iMC,:,:,:,:)./max(max(max(max(out(iMC,:,:,:,:))))), [min(min(min((rhoMC(:,:,:,iMC))))), max(max(max(rhoMC(:,:,:,iMC))))]);
        end
        cd(tmpPath);
        out = out .*conj(W); % .* (-2)
        clear 'rhoRef' 'sData';
    end
end
