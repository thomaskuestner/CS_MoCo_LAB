function dImg = ktFOCUSS3Dsingle(obj, kSpace)
% 3D FOCUSS algorithm (single precession) - OBSOLETE

% here: kSpace: k_y - k_z - k_x - cha
% calibSize: k_y - k_x - k_z
% fft: f -> t   x -> k_x   y -> k_y   z -> k_z

fprintf(' |- prepare k-space and image kernel...\n');

[nPha, nZ, nFreq, nCha] = size(kSpace);
% kSpace = fftshift(ifft(ifftshift(kSpace, 3),[],3),3).*nFreq; % -> k_y-k_z-x-cha space  ^= v (nü)
if(any(obj.trafo.fftBA(1,:)))
    for iBA=find(obj.trafo.fftBA(1,:))
        kSpace = ifftnshift(kSpace,iBA);
    end
end

% same subsampling for all channels
% lMask = repmat(lMask, [1 1 1 nCha]); % sparse matrix (CS pattern) => k_y - k_z - k_x/x - cha
% retrieve mask
% lMask = obj.fullMask; % k_y-k_z-k_x/x-cha

% calculate transformation basis
[fTrafo, bTrafo, kernelFTrafo, kernelBTrafo, paraTrafo] = compBasis(kSpace,obj);

% prepare initial estimate
fprintf(' |- prepare initial estimate...\n');
W = single(kSpace); % k_y-k_z-x-cha
helper = true(size(W));
m = [nPha, nZ, nFreq, nCha];
% find the k-space center
s = obj.calibSize; 
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
helper(idx{1},idx{2},idx{3},:) = false;
kSpaceCenter = ~helper;

% windowing for cut-out region
if(obj.window.windowingOn)
    window = windowND(obj.window.type,[length(idx{1}) length(idx{2})],obj.window.windowOpt{1},obj.window.windowOpt{2});
    window = repmat(window,[1 1 nFreq nCha]);
end

W(helper) = 0; % just take the k-space center (y values); all other values are zero
clear 'helper' 'idx' 'hShift' 'tmp' 'm' 's';
if(exist('window','var'))
    W(~helper) = W(~helper) .* window(:); 
    clear 'window';
end

pixEnergy = abs(W).^2;
totEnergy = squeeze(sum(sum(sum(pixEnergy,1),2),3));
W = abs(fTrafo(W)).^(obj.p); % rho values -> y-z-x-cha new_base
for c=1:nCha
%     helper = W(:,:,:,c);
%     W(:,:,:,c) = helper./max(helper(:)); % scaling
%     lScale = totEnergy(c) .* ones(size(W(:,:,:,c)));
%     tmp = repmat(mean(mean(pixEnergy(:,:,:,c),1),3),[size(lScale,1), 1, size(lScale,3)]);
%     lScale(pixEnergy(:,:,:,c) <= tmp) = 1./totEnergy(c);
%     W(:,:,:,c) = W(:,:,:,c)./lScale;
    W(:,:,:,c) = W(:,:,:,c)./totEnergy(c);
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
            kernelImg(:,:,:,:,iCha) = kernelFTrafo(fftnshift(obj.kernelImg(:,:,:,:,iCha),unique([obj.trafo.kspaceTrafoFFTdims, obj.trafo.fftdim]),unique([obj.trafo.kspaceTrafoScrambledims, obj.trafo.scrambledim])));
        else
            kernelImg(:,:,:,:,iCha) = kernelFTrafo(obj.kernelImg(:,:,:,:,iCha)); % => y-z-x-cha-cha (new_base)
        end
    end
    obj.kernelImg = []; % clear in object
end

fprintf(' |- start reconstruction...\n');

dispProgress('FOCUSS',0,obj.iNOUTER);
for iI = 1:obj.iNOUTER % outer loop for kt-FOCUSS
    % initial estimates for CG
    q = complex(zeros(size(W),'single'),zeros(size(W),'single'));
    rho = complex(zeros(size(W),'single'),zeros(size(W),'single'));
    d_old = complex(zeros(size(W),'single'),zeros(size(W),'single')); % with this intialization => d_0 = - g_0
    g_old = complex(ones (size(W),'single'),ones (size(W),'single')); % this initialization prevents division by zero
    dispProgress('CG',0,obj.iNINNER);
    for iJ = 1:obj.iNINNER % inner loop for conjugate gradient
        e = single(kSpace) - obj.fullMask.*bTrafo(rho); % ^= v - F*rho; kSpace (k_y-k_z-x-cha); rho (y-z-x-cha new_base); e (k_y-k_z-x-cha (cart))
        finished = false(1,nCha);
        for c=1:nCha
            helper = e(:,:,:,c);
            finished(c) = norm(helper(:)) <= obj.epsilon;
        end
        clear 'helper';
        if(all(finished))
            break;
        end
%         G = -conj(W).*fTrafo(e) .* nPha .* nZ + obj.lambda .*q + obj.lambdaCalib .* kernelMult(obj.lambdaCalib) + .5* obj.lambdaTV .* gradTV(obj.lambdaTV);
        G = -conj(W).*fTrafo(e) .* nPha .* nZ + obj.lambda .*q + obj.lambdaCalib .* kernelMult(obj.lambdaCalib) + .5* obj.lambdaTV .* gradTV(obj.lambdaTV) + obj.lambdaESPReSSo .* gradESPReSSO(obj.lambdaESPReSSo);
        beta = zeros(nPha,nZ,nFreq,nCha,'single');
        if(isa(G,'TRAFO'))
            beta = TRAFO(beta, getMeta(G), paraTrafo);
        end
        for c=1:nCha
            G_helper = G(:,:,:,c);
            g_old_helper = g_old(:,:,:,c);
%             beta(:,:,:,c) = (G_helper(:).'*G_helper(:))/(g_old_helper(:).'*g_old_helper(:)); % Fletcher-Reaves
            beta(:,:,:,c) = max(0,(G_helper(:).'*(G_helper(:)-g_old_helper(:)))/(g_old_helper(:).'*g_old_helper(:))); % Polak-Ribiere
%             d_old_helper = d_old(:,:,:,c); % Hestenes-Stiefel
%             beta(:,:,:,c) = -(G_helper(:).'*(G_helper(:)-g_old_helper(:)))/(d_old_helper(:).'*(G_helper(:)-g_old_helper(:)));
        end
        if(~isempty(paraTrafo.permRule)), beta = permute(beta,paraTrafo.permRule); end
        clear 'G_helper' 'g_old_helper';
        d = beta.*d_old - G;
        d_old = d;
        g_old = G;
        z = bTrafo(W.*d).*obj.fullMask; % -> k_y-k_z-x-cha
        alpha = zeros(nPha,nZ,nFreq,nCha,'single'); 
        if(isa(d,'TRAFO'))
            alpha = TRAFO(alpha, getMeta(d), paraTrafo);
        end
        for c=1:nCha
            z_helper = z(:,:,:,c);
            e_helper = e(:,:,:,c);
            alpha(:,:,:,c) = (z_helper(:)'*e_helper(:))/(z_helper(:)'*z_helper(:));
%             alpha(:,:,:,c) = abs((z_helper(:)'*e_helper(:))/(z_helper(:)'*z_helper(:))); % abs: to ensure a step into positive direction => wrong!
        end
        if(~isempty(paraTrafo.permRule)), alpha = permute(alpha,paraTrafo.permRule); end
        clear 'z_helper' 'e_helper';
        q = q + alpha.*d;
        rho = W.*q;
%         displayInterRes( kernelBTrafo(rho) );
        dispProgress('CG',iJ/obj.iNINNER);
        dispProgress('FOCUSS',((iI - 1)*obj.iNINNER + iJ)/(obj.iNOUTER*obj.iNINNER));
        if(obj.espresso.state && obj.espresso.reconType == 2)
    %     tmp = permute(abs(fft2shift(rho)) > 0,[1 3 2 4]);
            cs.smplPtrn = permute(obj.fullMask(:,:,:,1),[1 3 2]); %tmp(:,:,:,1);
            cs.pfn = obj.espresso.pfn;
            rho = bTrafo(rho);
            if(any(obj.trafo.fftBA(1,:))) % go back to kSpace
                for iBA=find(obj.trafo.fftBA(1,:))
                    rho = fftnshift(rho,iBA);
                end
            end
            rho = permute(espresso_recon(permute(rho,[4 1 3 2]), obj.espresso.iter, false, cs),[2 4 3 1]); % in: cha-k_y-k_x-k_z, out: y-z-x-cha (cart)
            if(obj.trafo.kspaceTrafo)
                rho = fftnshift(rho,unique([obj.trafo.kspaceTrafoFFTdims, obj.trafo.fftdim]),unique([obj.trafo.kspaceTrafoScrambledims, obj.trafo.scrambledim]));
            end
            rho = kernelFTrafo(rho); % y-z-x-cha (new_base)
        end
    end
    dispProgress('CG','Close');
    if(obj.espresso.state && obj.espresso.reconType == 2)
        cs.smplPtrn = permute(obj.fullMask(:,:,:,1),[1 3 2]); %tmp(:,:,:,1);
        cs.pfn = obj.espresso.pfn;
        rho = bTrafo(rho);
        if(any(obj.trafo.fftBA(1,:))) % go back to kSpace
            for iBA=find(obj.trafo.fftBA(1,:))
                rho = fftnshift(rho,iBA);
            end
        end
        rho = permute(espresso_recon(permute(rho,[4 1 3 2]), obj.espresso.iter, false, cs),[2 4 3 1]); % in: cha-k_y-k_x-k_z, out: y-z-x-cha (cart)
        if(obj.trafo.kspaceTrafo)
            rho = fftnshift(rho,unique([obj.trafo.kspaceTrafoFFTdims, obj.trafo.fftdim]),unique([obj.trafo.kspaceTrafoScrambledims, obj.trafo.scrambledim]));
        end
        rho = kernelFTrafo(rho); % y-z-x-cha (new_base)
    end
%     totEnergy = squeeze(sum(sum(sum(abs(bTrafo(rho)).^2,1),2),3));
    W = abs(rho).^(obj.p);
    for c=1:nCha
%         helper = W(:,:,:,c);
%         W(:,:,:,c) = helper./max(helper(:)); % scaling
        W(:,:,:,c) = W(:,:,:,c)./totEnergy(c);
    end
    clear 'helper';
end
dispProgress('FOCUSS', 'Close');
rho = W.*q;
dImg = kernelBTrafo(rho); % => y-z-x-cha (cart)
if(obj.trafo.kspaceTrafo)
    dImg = ifftnshift(dImg,obj.trafo.kspaceTrafoFFTdims,obj.trafo.kspaceTrafoScrambledims); % 3D: y and z
end
if(any(obj.trafo.fftBA(2,:)))
    for iBA=find(obj.trafo.fftBA(2,:))
        dImg = ifftnshift(dImg,iBA);
    end
end
dImg = permute(dImg, [1 3 2 4]); % => y-x-z-cha (cart)
dImg = shiftdim(mat2cell(dImg, nPha, nFreq, nZ, ones(1,nCha)),2);

    function res = kernelMult(lambdaCalib)
        if(lambdaCalib == 0)
            res = 0; % scalar needed for TRAFO object addition
            return;
        end
        res = zeros(size(kernelImg(:,:,:,:,1)));
        % zero-padding just allowed in cartesian image domain
        if(any(size(W) ~= size(res))) % same holds for q (q and W always have the same size)
            padsizeCurr = padsize(1:end-1);
            resize = true;
            Wtmp = kernelBTrafo(W);
            cropsize = size(Wtmp);
            if(obj.trafo.zeroPad)
                Wtmp = zpad(Wtmp,padsizeCurr);
            else
                Wout = zeros(padsizeCurr);
                RI = imref3d(size(Wtmp(:,:,:,1)),1,1,1);
                scaleFactor = size(Wout(:,:,:,1))./size(Wtmp(:,:,:,1));
                tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                for iCha=1:nCha
                    Wout(:,:,:,iCha) = crop(imwarp(Wtmp(:,:,:,iCha),RI,tFormResize,para.rescaleInterp),padsizeCurr(1:3));
                end
                Wtmp = Wout;
                clear 'Wout';
            end
            Wtmp = kernelFTrafo(Wtmp);

            qtmp = kernelBTrafo(q);
            if(obj.trafo.zeroPad)
                qtmp = zpad(qtmp,padsizeCurr);
            else
                qout = zeros(padsizeCurr);
                RI = imref3d(size(qtmp(:,:,:,1)),1,1,1);
                scaleFactor = size(qout(:,:,:,1))./size(qtmp(:,:,:,1));
                tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                for iCha=1:nCha
                    qout(:,:,:,iCha) = crop(imwarp(qtmp(:,:,:,iCha),RI,tFormResize,para.rescaleInterp),padsizeCurr(1:3));
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
            kernel = kernelImg(:,:,:,:,c); % y-z-x-cha
            kMI = (kernel - ones(size(kernel))).*Wtmp; %zpad(W,size(kernel)); % kernel - 1 % TODO: zpad in new basis incorrect > zpad in cart basis! -> btrafo - zpad - ftrafo
    %         res(:,:,:,c) = sum(kMI.*kMI .* q,4);
            res(:,:,:,c) = sum(conj(kMI) .* kMI .*qtmp,4); %zpad(q,size(kernel)),4);
        end
        clear 'Wtmp' 'qtmp' 'padsizeCurr' 'kMI' 'kernel';

        % res = crop(res,size(q));
        if(resize)
            res = kernelBTrafo(res);
            if(obj.trafo.zeroPad)
                res = crop(res,cropsize); 
            else
                resout = zeros(cropsize);
                RI = imref3d(size(res(:,:,:,1)),1,1,1);
                scaleFactor = size(resout(:,:,:,1))./size(res(:,:,:,1));
                tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                for iCha=1:nCha
                    resout(:,:,:,iCha) = crop(imwarp(res(:,:,:,iCha),RI,tFormResize,para.rescaleInterp),cropsize(1:3));
                end
                res = resout;
                clear 'resout';
            end
            res = kernelFTrafo(res);
        end
    end

    function tv = gradTV(lambdaTV)
        % calculate gradient of TV from image
        if(lambdaTV == 0)
            tv = 0;
            return;
        end
        rhotmp = kernelBTrafo(W.*q);
        if(nnz(rhotmp) == 0)
            tv = 0; % scalar needed for TRAFO object addition
            return;
        end
        tv = zeros(size(rhotmp));
        type = 2;
        l1Smooth = 1e-15;
        p = 1;

        for c=1:nCha
            x = rhotmp(:,:,:,c);     
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

                tv(:,:,:,c) = adjDx + adjDy + adjDz;
                

            elseif(type == 2)
                % Matlab calculated differentiation
                % corner points
                i = 1; j = 1; k = 1;
                tv(i,j,k,c) = -(x(i + 1, j, k) - 3*x(i, j, k) + x(i, j, k + 1) + x(i, j + 1, k))./((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2);
                i = size(tv,1); j = 1; k = 1;
                tv(i,j,k,c) = (x(i, j, k) - x(i - 1, j, k))./((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2) - (2*x(i, j, k + 1) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2).^(1./2));
                i = 1; j = size(tv,2); k = 1;
                tv(i,j,k,c) = (x(i, j, k) - x(i, j - 1, k))./((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                i = size(tv,1); j = size(tv,2); k = 1;
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2).^(1./2));
                i = 1; j = 1; k = size(tv,3);
                tv(i,j,k,c) = (x(i, j, k) - x(i, j, k - 1))./((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                i = size(tv,1); j = 1; k = size(tv,3);
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2).^(1./2));
                i = 1; j = size(tv,2); k = size(tv,3);
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i + 1, j, k))./(2*((x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2));
                i = size(tv,1); j = size(tv,2); k = size(tv,3);
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j, k)).^2).^(1./2));
                % outer lanes
                i = 2:size(tv,1)-1; j = 1; k = 1;
                tv(i,j,k,c) = (x(i, j, k) - x(i - 1, j, k))./((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2) - (x(i + 1, j, k) - 3*x(i, j, k) + x(i, j, k + 1) + x(i, j + 1, k))./((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2);
                i = 1; j = 2:size(tv,2)-1; k = 1;
                tv(i,j,k,c) = (x(i, j, k) - x(i, j - 1, k))./((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2) - (x(i + 1, j, k) - 3*x(i, j, k) + x(i, j, k + 1) + x(i, j + 1, k))./((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2);
                i = size(tv,1); j = 2:size(tv,2)-1; k = 1;
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (x(i, j, k) - x(i - 1, j, k))./((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2) - (2*x(i, j, k + 1) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2).^(1./2));
                i = 2:size(tv,1)-1; j = size(tv,2); k = 1;
                tv(i,j,k,c) = (x(i, j, k) - x(i, j - 1, k))./((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                i = 2:size(tv,1)-1; j = 1; k = size(tv,3);
                tv(i,j,k,c) = (x(i, j, k) - x(i, j, k - 1))./((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                i = 2:size(tv,1)-1; j = size(tv,2); k = size(tv,3);
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i + 1, j, k))./(2*((x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2));
                i = 1; j = 2:size(tv,2)-1; k = size(tv,3);
                tv(i,j,k,c) = (x(i, j, k) - x(i, j, k - 1))./((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2));
                i = size(tv,1); j = 2:size(tv,2)-1; k = size(tv,3);
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2).^(1./2));
                i = 1; j = 1; k = 2:size(tv,3)-1;
                tv(i,j,k,c) = (x(i, j, k) - x(i, j, k - 1))./((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2) - (x(i + 1, j, k) - 3*x(i, j, k) + x(i, j, k + 1) + x(i, j + 1, k))./((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2);
                i = 1; j = size(tv,2); k = 2:size(tv,3)-1;
                tv(i,j,k,c) = (x(i, j, k) - x(i, j - 1, k))./((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2));
                i = size(tv,1); j = 1; k = 2:size(tv,3)-1;
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (x(i, j, k) - x(i - 1, j, k))./((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2) - (2*x(i, j, k + 1) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2).^(1./2));
                i = size(tv,1); j = size(tv,2); k = 2:size(tv,3)-1;
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2).^(1./2));
                % outer planes
                i = 2:size(tv,1)-1; j = 2:size(tv,2)-1; k = 1;
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 6*x(i, j, k) + 2*x(i, j, k + 1) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                i = 2:size(tv,1)-1; j = 2:size(tv,2)-1; k = size(tv,3);
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2));
                i = 1; j = 2:size(tv,2)-1; k = 2:size(tv,3)-1;
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 6*x(i, j, k) + 2*x(i, j, k + 1) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                i = size(tv,1); j = 2:size(tv,2)-1; k = 2:size(tv,3)-1;
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i, j, k + 1) - 4*x(i, j, k) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2));
                i = 2:size(tv,1)-1; j = 1; k = 2:size(tv,3)-1;
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 6*x(i, j, k) + 2*x(i, j, k + 1) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));
                i = 2:size(tv,1)-1; j = size(tv,2); k = 2:size(tv,3)-1;
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 4*x(i, j, k) + 2*x(i, j, k + 1))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2));
                % inner points
                i = 2:size(tv,1)-1; j = 2:size(tv,2)-1; k = 2:size(tv,3)-1;
                tv(i,j,k,c) = (2*x(i, j, k) - 2*x(i, j, k - 1))./(2*((x(i, j, k - 1) - x(i, j + 1, k - 1)).^2 + (x(i, j, k - 1) - x(i + 1, j, k - 1)).^2 + (x(i, j, k - 1) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i, j - 1, k))./(2*((x(i, j - 1, k) - x(i, j - 1, k + 1)).^2 + (x(i, j - 1, k) - x(i + 1, j - 1, k)).^2 + (x(i, j - 1, k) - x(i, j, k)).^2).^(1./2)) + (2*x(i, j, k) - 2*x(i - 1, j, k))./(2*((x(i - 1, j, k) - x(i - 1, j, k + 1)).^2 + (x(i - 1, j, k) - x(i - 1, j + 1, k)).^2 + (x(i - 1, j, k) - x(i, j, k)).^2).^(1./2)) - (2*x(i + 1, j, k) - 6*x(i, j, k) + 2*x(i, j, k + 1) + 2*x(i, j + 1, k))./(2*((x(i, j, k + 1) - x(i, j, k)).^2 + (x(i, j + 1, k) - x(i, j, k)).^2 + (x(i + 1, j, k) - x(i, j, k)).^2).^(1./2));





            else
    %         % corner points
    %         tv(1,1,1,c) = (-(x(2,1,1) - x(1,1,1)) - (x(1,2,1) - x(1,1,1)) - (x(1,1,2) - x(1,1,1)))/sqrt((x(2,1,1) - x(1,1,1)).^2 + (x(1,2,1) - x(1,1,1)).^2 + (x(1,1,2) - x(1,1,1)).^2);
    %         tv(end,1,1,c) = (x(end,1,1) - x(end-1,1,1))/sqrt((x(end,1,1) - x(end-1,1,1)).^2 + (x(end-1,1,2) - x(end-1,1,2)).^2 + (x(end-1,2,1) - x(end-1,1,1)).^2) - ((x(end,1,2) - 2*x(end,1,1) + x(end,2,1))/sqrt((x(end,1,2) - x(end,1,1)).^2 + (x(end,2,1) - x(end,1,1)).^2));
    %         tv(1,end,1,c) = (-(x(1,end-1,1) - x(1,end,1)))/sqrt((x(1,end-1,1) - x(1,end-1,2)).^2 + (x(1,end-1,1) - x(2,end-1,1)).^2 + (x(1,end-1,1) - x(1,end,1)).^2) - ((x(1,end,2) - 2*x(1,end,1) + x(2,end,1))/sqrt((x(1,end,1) - x(1,end,2)).^2 + (x(1,end,1) - x(2,end,1)).^2));
    %         tv(end,end,1,c) = (-(x(end-1,end,1) - x(end,end,1)))/sqrt((x(end-1,end,1) - x(end-1,end,2)).^2 + (x(end-1,end,1) - x(end,end,1)).^2) - ((x(end,end-1,1) - x(end,end,1))/sqrt((x(end,end-1,1) - x(end,end-1,2)).^2 + (x(end,end-1,1) - x(end,end,1)).^2));
    %         tv(1,1,end,c) = (-(x(1,1,end-1) - x(1,1,end)))/sqrt((x(1,1,end-1) - x(1,1,end)).^2 + (x(1,1,end-1) - x(1,2,end-1)).^2 + (x(1,1,end-1) - x(2,1,end-1)).^2) - ((x(1,2,end) - 2*x(1,1,end) + x(2,1,end))/sqrt((x(1,1,end) - x(1,2,end)).^2 + (x(1,1,end) - x(2,1,end)).^2));
    %         tv(end,1,end,c) = (-(x(end,1,end-1) - x(end,1,end)))/sqrt((x(end,1,end-1) - x(end,1,end)).^2 + (x(end,1,end-1) - x(end,2,end-1)).^2) - ((x(end-1,1,end) - x(end,1,end))/sqrt((x(end-1,1,end) - x(end-1,2,end)).^2 + (x(end-1,1,end) - x(end,1,end)).^2));
    %         tv(1,end,end,c) = (-(x(1,end,end-1) - x(1,end,end)))/sqrt((x(1,end,end-1) - x(1,end,end)).^2 + (x(1,end,end-1) - x(2,end,end-1)).^2) - ((x(1,end-1,end) - x(1,end,end))/sqrt((x(1,end-1,end) - x(1,end,end)).^2 + (x(1,end-1,end) - x(2,end-1,end)).^2));
    %         tv(end,end,end,c) = 0;
    %         
    %         % outer lanes
    %         i = 2:size(tv,1)-1; j = 1; k = 1;
    %         tv(i,j,k,c) = (x(i,j,k) - x(i-1,j,k))./sqrt((x(i,j,k) - x(i-1,j,k)).^2 + (x(i-1,j+1,k) - x(i-1,j,k)).^2 + (x(i-1,j,k+1) - x(i-1,j,k)).^2) + (-(x(i+1,j,k) - x(i,j,k)) - (x(i,j+1,k) - x(i,j,k)) - (x(i,j,k+1) - x(i,j,k)))./sqrt((x(i+1,j,k) - x(i,j,k)).^2 + (x(i,j+1,k) - x(i,j,k)).^2 + (x(i,j,k+1) - x(i,j,k)).^2);
    %         i = 1; j = 2:size(tv,2)-2; k = 1;
    %         tv(i,j,k,c) = (-(x(i,j-1,k) - x(i,j,k)))./sqrt((x(i,j-1,k) - x(i,j-1,k+1)).^2 + (x(i,j-1,k) - x(i,j,k)).^2 + (x(i,j-1,k) - x(i+1,j-1,k)).^2) - (x(i,j,k+1) - 3*x(i,j,k) + x(i,j+1,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i,j+1,k)).^2 + (x(i,j,k) - x(i+1,j,k)).^2);
    %         i = size(tv,1); j = 2:size(tv,2)-1; k = 1;
    %         tv(i,j,k,c) = (-(x(i-1,j,k) - x(i,j,k)))./sqrt((x(i-1,j,k) - x(i-1,j,k+1)).^2 + (x(i-1,j,k) - x(i-1,j+1,k)).^2 + (x(i-1,j,k) - x(i,j,k)).^2) - (x(i,j,k+1) - 2*x(i,j,k) + x(i,j+1,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i,j+1,k)).^2) - (x(i,j-1,k) - x(i,j,k))./sqrt((x(i,j-1,k) - x(i,j-1,k+1)).^2 + (x(i,j-1,k) - x(i,j,k)).^2);
    %         i = 2:size(tv,1)-1; j = size(tv,2); k = 1;
    %         tv(i,j,k,c) = (-(x(i,j-1,k) - x(i,j,k)))./sqrt((x(i,j-1,k) - x(i,j-1,k+1)).^2 + (x(i,j-1,k) - x(i+1,j-1,k)).^2 + (x(i,j-1,k) - x(i,j,k)).^2) - (x(i,j,k+1) - 2*x(i,j,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i+1,j,k)).^2) - (x(i-1,j,k) - x(i,j,k))./sqrt((x(i-1,j,k) - x(i-1,j,k+1)).^2 + (x(i-1,j,k) - x(i,j,k)).^2);
    %         i = 2:size(tv,1)-1; j = 1; k = size(tv,3);
    %         tv(i,j,k,c) = (-(x(i,j,k-1) - x(i,j,k)))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i,j+1,k-1)).^2 + (x(i,j,k-1) - x(i+1,j,k-1)).^2) - (x(i,j+1,k) - 2*x(i,j,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j+1,k)).^2 + (x(i,j,k) - x(i+1,j,k)).^2) - (x(i-1,j,k) - x(i,j,k))./sqrt((x(i-1,j,k) - x(i-1,j+1,k)).^2 + (x(i-1,j,k) - x(i,j,k)).^2);
    %         i = 2:size(tv,1)-1; j = size(tv,2); k = size(tv,3);
    %         tv(i,j,k,c) = (-(x(i,j,k-1) - x(i,j,k)))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i+1,j,k-1)).^2) - (x(i,j-1,k) - x(i,j,k))./sqrt((x(i,j-1,k) - x(i+1,j-1,k)).^2 + (x(i,j-1,k) - x(i,j,k)).^2);
    %         i = 1; j = 2:size(tv,2)-1; k = size(tv,3);
    %         tv(i,j,k,c) = (-(x(i,j,k-1) - x(i,j,k)))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i,j+1,k-1)).^2 + (x(i,j,k-1) - x(i+1,j,k-1)).^2) - (x(i,j+1,k) - 2*x(i,j,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j+1,k)).^2 + (x(i,j,k) - x(i+1,j,k)).^2) - (x(i,j-1,k) - x(i,j,k))./sqrt((x(i,j-1,k) - x(i,j,k)).^2 + (x(i,j-1,k) - x(i+1,j-1,k)).^2);
    %         i = size(tv,1); j = 2:size(tv,2)-1; k = size(tv,3);
    %         tv(i,j,k,c) = (-(x(i,j,k-1) - x(i,j,k)))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i,j+1,k-1)).^2) - (x(i-1,j,k) - x(i,j,k))./sqrt((x(i-1,j,k) - x(i-1,j+1,k)).^2 + (x(i-1,j,k) - x(i,j,k)).^2);
    %         i = 1; j = 1; k = 2:size(tv,3)-1;
    %         tv(i,j,k,c) = (-(x(i,j,k-1) - x(i,j,k)))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i,j+1,k-1)).^2 + (x(i,j,k-1) - x(i+1,j,k-1)).^2) - (x(i,j,k+1) - 3*x(i,j,k) + x(i,j+1,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i,j+1,k)).^2 + (x(i,j,k) - x(i+1,j,k)).^2);
    %         i = 1; j = size(tv,2); k = 2:size(tv,3)-1;
    %         tv(i,j,k,c) = (-(x(i,j-1,k) - x(i,j,k)))./sqrt((x(i,j-1,k) - x(i,j-1,k+1)).^2 + (x(i,j-1,k) - x(i+1,j-1,k)).^2 + (x(i,j-1,k) - x(i,j,k)).^2) - (x(i,j,k+1) - 2*x(i,j,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i+1,j,k)).^2) - (x(i,j,k-1) - x(i,j,k))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i+1,j,k-1)).^2);
    %         i = size(tv,1); j = 1; k = 2:size(tv,3)-1;
    %         tv(i,j,k,c) = (-(x(i-1,j,k) - x(i,j,k)))./sqrt((x(i-1,j,k) - x(i-1,j,k+1)).^2 + (x(i-1,j,k) - x(i-1,j+1,k)).^2 + (x(i-1,j,k) - x(i,j,k)).^2) - (x(i,j,k+1) - 2*x(i,j,k) + x(i,j+1,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i,j+1,k)).^2) - (x(i,j,k-1) - x(i,j,k))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i,j+1,k-1)).^2);
    %         i =size(tv,1); j = size(tv,2); k = 2:size(tv,3)-1;
    %         tv(i,j,k,c) = (-(x(i-1,j,k) - x(i,j,k)))./sqrt((x(i-1,j,k) - x(i-1,j,k+1)).^2 + (x(i-1,j,k) - x(i,j,k)).^2) - (x(i,j-1,k) - x(i,j,k))./sqrt((x(i,j-1,k) - x(i,j-1,k+1)).^2 + (x(i,j-1,k) - x(i,j,k)).^2);
    %         
    %         % inner points
    %         i = 2:size(tv,1)-1; j = 2:size(tv,2)-1; k = 2:size(tv,3)-1;
    %         tv(i,j,k,c) = (-(x(i-1,j,k) - x(i,j,k)))./sqrt((x(i-1,j,k) - x(i-1,j,k+1)).^2 + (x(i-1,j,k) - x(i-1,j+1,k)).^2 + (x(i-1,j,k) - x(i,j,k)).^2) - (x(i,j-1,k) - x(i,j,k))./sqrt((x(i,j-1,k) - x(i,j-1,k+1)).^2 + (x(i,j-1,k) - x(i,j,k)).^2 + (x(i,j-1,k) - x(i+1,j-1,k)).^2) - (x(i,j,k-1) - x(i,j,k))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i,j+1,k-1)).^2 + (x(i,j,k-1) - x(i+1,j,k-1)).^2) - (x(i,j,k+1) - 3*x(i,j,k) + x(i,j+1,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i,j+1,k)).^2 + (x(i,j,k) - x(i+1,j,k)).^2);
    %         
    %         % outer planes
    %         i = 2:size(tv,1)-1; j = 2:size(tv,2)-1; k = 1;
    %         tv(i,j,k,c) = (-(x(i-1,j,k) - x(i,j,k)))./sqrt((x(i-1,j,k) - x(i-1,j,k+1)).^2 + (x(i-1,j,k) - x(i-1,j+1,k)).^2 + (x(i-1,j,k) - x(i,j,k)).^2) - (x(i,j-1,k) - x(i,j,k))./sqrt((x(i,j-1,k) - x(i,j-1,k+1)).^2 + (x(i,j-1,k) - x(i,j,k)).^2 + (x(i,j-1,k) - x(i+1,j-1,k)).^2) - (x(i,j,k+1) - 3*x(i,j,k) + x(i,j+1,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i,j+1,k)).^2 + (x(i,j,k) - x(i+1,j,k)).^2);
    %         i = 2:size(tv,1)-1; j = 2:size(tv,2)-1; k = size(tv,3);
    %         tv(i,j,k,c) = (-(x(i,j,k-1) - x(i,j,k)))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i,j+1,k-1)).^2 + (x(i,j,k-1) - x(i+1,j,k-1)).^2) - (x(i,j+1,k) - 2*x(i,j,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j+1,k)).^2 + (x(i,j,k) - x(i+1,j,k)).^2) - (x(i-1,j,k) - x(i,j,k))./sqrt((x(i-1,j,k) - x(i-1,j+1,k)).^2 + (x(i-1,j,k) - x(i,j,k)).^2) - (x(i,j-1,k) - x(i,j,k))./sqrt((x(i,j-1,k) - x(i,j,k)).^2 + (x(i,j-1,k) - x(i+1,j-1,k)).^2);
    %         i = 1; j = 2:size(tv,2)-1; k = 2:size(tv,3)-1;
    %         tv(i,j,k,c) = (-(x(i,j-1,k) - x(i,j,k)))./sqrt((x(i,j-1,k) - x(i,j-1,k+1)).^2 + (x(i,j-1,k) - x(i,j,k)).^2 + (x(i,j-1,k) - x(i+1,j-1,k)).^2) - (x(i,j,k-1) - x(i,j,k))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i,j+1,k-1)).^2 + (x(i,j,k-1) - x(i+1,j,k-1)).^2) - (x(i,j,k+1) - 3*x(i,j,k) + x(i,j+1,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i,j+1,k)).^2 + (x(i,j,k) - x(i+1,j,k)).^2);
    %         i = size(tv,1); j = 2:size(tv,2)-1; k = 2:size(tv,3)-1;
    %         tv(i,j,k,c) = (-x(i-1,j,k) - x(i,j,k))./sqrt((x(i-1,j,k) - x(i-1,j,k+1)).^2 + (x(i-1,j,k) - x(i-1,j+1,k)).^2 + (x(i-1,j,k) - x(i,j,k)).^2) - (x(i,j,k+1) - 2*x(i,j,k) + x(i,j+1,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i,j+1,k)).^2) - (x(i,j-1,k) - x(i,j,k))./sqrt((x(i,j-1,k) - x(i,j-1,k+1)).^2 + (x(i,j-1,k) - x(i,j,k)).^2) - (x(i,j,k-1) - x(i,j,k))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i,j+1,k-1)).^2);
    %         i = 2:size(tv,1)-1; j = 1; k = 2:size(tv,3)-1;
    %         tv(i,j,k,c) = (-(x(i-1,j,k) - x(i,j,k)))./sqrt((x(i-1,j,k) - x(i-1,j,k+1)).^2 + (x(i-1,j,k) - x(i-1,j+1,k)).^2 + (x(i-1,j,k) - x(i,j,k)).^2) - (x(i,j,k-1) - x(i,j,k))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i,j+1,k-1)).^2 + (x(i,j,k-1) - x(i+1,j,k-1)).^2) - (x(i,j,k+1) - 3*x(i,j,k) + x(i,j+1,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i,j+1,k)).^2 + (x(i,j,k) - x(i+1,j,k)).^2);
    %         i = 2:size(tv,1)-1; j = size(tv,2); k = 2:size(tv,3)-1;
    %         tv(i,j,k,c) = (-(x(i,j-1,k) - x(i,j,k)))./sqrt((x(i,j-1,k) - x(i,j-1,k+1)).^2 + (x(i,j-1,k) - x(i+1,j-1,k)).^2 + (x(i,j-1,k) - x(i,j,k)).^2) - (x(i,j,k+1) - 2*x(i,j,k) + x(i+1,j,k))./sqrt((x(i,j,k) - x(i,j,k+1)).^2 + (x(i,j,k) - x(i+1,j,k)).^2) - (x(i-1,j,k) - x(i,j,k))./sqrt((x(i-1,j,k) - x(i-1,j,k+1)).^2 + (x(i-1,j,k) - x(i,j,k)).^2) - (x(i,j,k-1) - x(i,j,k))./sqrt((x(i,j,k-1) - x(i,j,k)).^2 + (x(i,j,k-1) - x(i+1,j,k-1)).^2);        
            end
        end
        tv = kernelFTrafo(tv);
        clear 'Gtv' 'adjDx' 'adjDy' 'adjDz' 'tmp';
    end

    function out = gradESPReSSO(lambdaESPReSSo)
        % calculate ESPReSSo gradient
        if(lambdaESPReSSo == 0 || obj.espresso.reconType ~= 1)
            out = 0;
            return;
        end
        
        if(~obj.espresso.state)
%             mask = obj.fullMask(:,:,:,1); % y-z-x
            maskRight = obj.fullMask;
            maskLeft = maskRight;
            maskRight(1:end/2,:,:,:) = false;
            maskLeft(end/2:end,:,:,:) = false;
            maskRight = maskRight | kSpaceCenter; % y-z-x-cha
            maskLeft = maskLeft | kSpaceCenter;
            cs.smplPtrn = permute(maskRight(:,:,:,1),[1 3 2]);
            cs.pfn = 0.5;
        else      
            cs.smplPtrn = permute(obj.fullMask(:,:,:,1),[1 3 2]); %tmp(:,:,:,1);
            cs.pfn = obj.espresso.pfn;
        end
        
        rhoTmp = bTrafo(W.*q);
        if(nnz(rhoTmp) == 0)
            rhoTmp = kSpace;
        end
        if(any(obj.trafo.fftBA(1,:))) % go back to kSpace
            for iBA=find(obj.trafo.fftBA(1,:))
                rhoTmp = fftnshift(rhoTmp,iBA);
            end
        end
        if(nnz(rhoTmp) > 0)
            [lmaskSym,pfDim,phase] = espresso_reconGrad(permute(rhoTmp,[4 1 3 2]), cs); % in: cha-k_y-k_x-k_z, out: y-z-x-cha (cart)
        else
            [lmaskSym,pfDim,phase] = espresso_reconGrad(permute(rhoTmp,[4 1 3 2]), cs); % in: cha-k_y-k_x-k_z, out: y-z-x-cha (cart)
            phase = zeros(size(phase));
        end

%         if(obj.trafo.kspaceTrafo)
%             rhoTmp = fftnshift(rhoTmp,unique([obj.trafo.kspaceTrafoFFTdims, obj.trafo.fftdim]),unique([obj.trafo.kspaceTrafoScrambledims, obj.trafo.scrambledim]));
%         end
        phase = permute(phase,[1 3 2 4]);
        phase = exp(2i * angle(phase));
        lmaskSym = repmat(squeeze(lmaskSym),[1 1 1 nCha]);
        lmaskSym = permute(lmaskSym,[1 3 2 4]);
        if(~obj.espresso.state)
            lmaskConj = maskRight(end:-1:1,:,:,:);
            
            kSpaceCombi = kSpace;
            kSpaceCombi(lmaskConj) = conj(kSpace(maskRight));
            lmaskConj = maskLeft(end:-1:1,:,:);
            kSpaceCombi(lmaskConj) = conj(kSpace(maskLeft));
            
        else
            lmaskConj = xor(lmaskSym,obj.fullMask);
        
            if(pfDim == 1) %y 
                lmaskLower = lmaskConj(end:-1:1,:,:,:);
            elseif(pfDim == 2) % x
                lmaskLower = lmaskConj(:,:,end:-1:1,:);
            elseif(pfDim == 3) % z
                lmaskLower = lmaskConj(:,end:-1:1,:,:);
            end
        
            kSpaceCombi = kSpace;
            kSpaceCombi(lmaskConj) = conj(kSpace(lmaskLower));
        end

        out = 0.5 * (conj(W) .* q .* (1+phase) + conj(W) .* W .* conj(q) .* phase - fTrafo(kSpaceCombi) .* (1+phase));        
    end
end