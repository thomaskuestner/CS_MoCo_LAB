function dImg = CSNLCG(obj, kSpace)
% non-linear conjugate gradient

% here: kSpace: k_y - k_z - k_x - cha

% ||v - Phi*F*rho||_2^2 + lambda * ||Psi * rho||_1 
% ATTENTION: directions of Psi are changed compared to FOCUSS
% definition!!!
% rho is in cartesian basis and not sparse
% Psi * rho is in sparse basis
%
% input:
% obj           CS reconstruction object (holding all parameters)
% kSpace        subsampled kSpace
%
% output:
% dImg          reconstructed channel individual image
%
% (c) Thomas Kuestner
% -------------------------------------------------------------------------

fprintf(' |- prepare k-space and image kernel...\n');

[nPha, nZ, nFreq, nCha] = size(kSpace);
% kSpace = fftshift(ifft(ifftshift(kSpace, 3),[],3),3).*nFreq; % -> k_y-k_z-x-cha space  ^= v (nü)
if(any(obj.trafo.fftBA(1,:)))
    for iBA=find(obj.trafo.fftBA(1,:))
        kSpace = ifftnshift(kSpace,iBA);
    end
end

% prepare initial estimate
fprintf(' |- prepare initial estimate...\n');
W = kSpace; % k_y-k_z-x-cha
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
if(length(idx) == 2), idx{3} = 1; end;
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

%%% OR: take 
% rho = fTrafo(kSpace);

if(obj.lSense)
    dImg = W; % k_y-k_z-x-cha
    iTrafodim = find(obj.trafo.fftBA(1,:) == false);
    iTrafodim = iTrafodim(ismember(iTrafodim,[1 2 3])); % except for t
    for iBA=iTrafodim
        dImg = ifftnshift(dImg,iBA);
    end
    clear 'iTrafodim' 'iBA';
    dImgCenter = dImg;
    dImg = sqrt(sum(abs(dImg).^2,4)); 
    dSensemap = dImgCenter./(repmat(dImg,[1 1 1 size(dImgCenter,4)])); % y-x-z-cha
    clear 'dImgCenter' 'dImg';
    dWindow = windowND(@hamming, [size(dSensemap,1), size(dSensemap,2), size(dSensemap,3)]);
    dSensemap = fftnshift(dSensemap,1:3);
    dSensemap = dSensemap .* repmat(dWindow, [1 1 1 size(dSensemap,4) size(dSensemap,5)]);
    dSensemap = ifftnshift(dSensemap,1:3);
    dSensemapFctr = sqrt(sum(conj(dSensemap).*dSensemap,4)).^(-1);
    dSensemapFctr(isinf(dSensemapFctr)) = 0;
    dSensemap = conj(repmat(dSensemapFctr,[1 1 1 size(dSensemap,4)]) .* dSensemap);
    clear 'dSensemapFctr' 'dWindow';
else
    dSensemap = 1;
end

% calculate transformation basis
[fTrafo, bTrafo, kernelFTrafo, kernelBTrafo, paraTrafo] = compBasis(kSpace,obj,dSensemap);

% rho = fTrafo(W);
rho = complex(zeros(size(W),obj.measPara.precision),zeros(size(W),obj.measPara.precision));

% line search parameters
maxlsiter = obj.lineSearchItnlim ;
gradToll = obj.gradToll ;
alpha = obj.lineSearchAlpha; 
beta = obj.lineSearchBeta;
t0 = obj.lineSearchT0;
k = 0;

g_old = objgrad(rho);
dx = -g_old;

% if(obj.lambdaCalib > 0)
%     % prepare image kernel
%     kernelImg = zeros(size(obj.kernelImg),obj.measPara.precision);
%     padsize = size(kernelImg);
%     if(isa(W,'TRAFO'))
%         kernelImg = TRAFO(kernelImg, getMeta(W), paraTrafo, 'kernel');
%     end
%     for iCha=1:nCha
%         if(obj.trafo.kspaceTrafo)           
%             kernelImg(:,:,:,:,iCha) = kernelFTrafo(fftnshift(obj.kernelImg(:,:,:,:,iCha),unique([obj.trafo.kspaceTrafoFFTdims, obj.trafo.fftdim]),unique([obj.trafo.kspaceTrafoScrambledims, obj.trafo.scrambledim])));
%         else
%             kernelImg(:,:,:,:,iCha) = kernelFTrafo(obj.kernelImg(:,:,:,:,iCha)); % => y-z-x-cha-cha (new_base)
%         end
%     end
%     obj.kernelImg = []; % clear in object
% end

fprintf(' |- start reconstruction...\n');

dispProgress('CG',0,obj.iNINNER);

while((k < obj.iNINNER) && (norm(dx(:)) > gradToll)) % non-linear CG
    % backtracking line-search
    f0 = objective(rho,dx,0);
	t = t0;
    f1 = objective(rho,dx,t);
	lsiter = 0;
    dispProgress('Line Search',0,maxlsiter)
    while (f1 > f0 - alpha*t*abs(g_old(:)'*dx(:)).^2 && (lsiter<maxlsiter))
		lsiter = lsiter + 1;
		t = t * beta;
		f1 = objective(rho,dx,t);
        dispProgress('Line Search',lsiter/maxlsiter);
    end
    dispProgress('Line Search','close');
    
    if lsiter == maxlsiter
		disp('Error - line search ...');
		return;
    end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2, t0 = t0 * beta;end 
	if lsiter<1, t0 = t0 / beta; end

    % update rho
	rho = (rho + t*dx);
    
    %conjugate gradient calculation
	G = objgrad(rho);
	bk = G(:)'*G(:)/(g_old(:)'*g_old(:)+eps);
	g_old = G;
	dx =  - G + bk* dx;
	k = k + 1;
    dispProgress('CG',k/obj.iNINNER);
end
dispProgress('CG','close');
  
dImg = rho; % => y-z-x-cha (cart)
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

    function res = objgrad(rho)
        % L2 norm part
        L2Grad = 2.*(ifftnshift(kSpace - obj.fullMask.*fftnshift(rho,obj.trafo.fftdim,obj.trafo.scrambledim),obj.trafo.fftdim,obj.trafo.scrambledim));
        
        % L1 norm part
        % sparsify
        w = kernelFTrafo(rho);
        Gsparse = kernelBTrafo(w.* (w.*conj(w)+obj.l1Smooth).^(-0.5));
        
        % TV
        if(obj.lambdaTV)
            Dx = rho([2:end,end],:,:,:) - rho;
            Dy = rho(:,[2:end,end],:,:) - rho;
            Dz = rho(:,:,[2:end,end],:) - rho;

            res = cat(5,Dx,Dy,Dz);
            Gtv = res.*(res.*conj(res) + obj.l1Smooth).^(-0.5);
            
            tmp = Gtv(:,:,:,:,1);
            adjDx = tmp([1,1:end-1],:,:,:) - tmp;
            adjDx(1,:,:,:) = -tmp(1,:,:,:);
            adjDx(end,:,:,:) = tmp(end-1,:,:,:);
            tmp = Gtv(:,:,:,:,2);
            adjDy = tmp(:,[1,1:end-1],:,:) - tmp;
            adjDy(:,1,:,:) = -tmp(:,1,:,:);
            adjDy(:,end,:,:) = tmp(:,end-1,:,:);
            tmp = Gtv(:,:,:,:,3);
            if(size(tmp,3) > 1)
                adjDz = tmp(:,:,[1,1:end-1],:) - tmp;
                adjDz(:,:,1,:) = -tmp(:,:,1,:);
                adjDz(:,:,end,:) = tmp(:,:,end-1,:);
            else
                adjDz = 0;
            end
            Gtv = adjDx + adjDy + adjDz;
        else
            Gtv = 0;
        end
        
        % ESPReSSo
        if(obj.lambdaESPReSSo)
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
        
            if(obj.espresso.constraint == 1) % POCS
                rhoTmp = fftnshift(rho,obj.trafo.fftdim,obj.trafo.scrambledim);
                if(any(obj.trafo.fftBA(1,:))) % go back to kSpace
                    for iBA=find(obj.trafo.fftBA(1,:))
                        rhoTmp = fftnshift(rhoTmp,iBA);
                    end
                end
                [lmaskSym,pfDim,phase] = espresso_reconGrad(permute(rhoTmp,[4 1 3 2]), cs); % in: cha-k_y-k_x-k_z, out: y-z-x-cha (cart)
                phase = permute(phase,[1 3 2 4]);
                phase = exp(2i * angle(phase));
                lmaskSym = repmat(squeeze(lmaskSym),[1 1 1 nCha]);
%                 lmaskSym = permute(lmaskSym,[1 3 2 4]);
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
                
                if(obj.espresso.norm == 1) % L1
                    l1Smooth = 1e-15;
                    rhoTmp = fftnshift(conj(rho).*phase + rho,obj.trafo.fftdim,obj.trafo.scrambledim);
                    if(any(obj.trafo.fftBA(1,:))) % go back to kSpace
                        for iBA=find(obj.trafo.fftBA(1,:))
                            rhoTmp = fftnshift(rhoTmp,iBA);
                        end
                    end
                    Gesp = rhoTmp - kSpaceCombi;
                    Gesp = Gesp(:)'*conj(Gesp(:));
                	Gesp = (Gesp+l1Smooth).^(-0.5) .* ( 2 .* conj(rho) + rho .* conj(phase) - ifftnshift(kSpaceCombi,obj.trafo.fftdim,obj.trafo.scrambledim) .* conj(phase));
                    % or: Gesp = ifftnshift(real(Gesp).* (Gesp.*conj(Gesp)+obj.l1Smooth).^(-0.5),obj.trafo.fftdim,obj.trafo.scrambledim);
                else % L2       
                    Gesp = 2 .* conj(rho) + rho .* conj(phase) - ifftnshift(kSpaceCombi,obj.trafo.fftdim,obj.trafo.scrambledim) .* conj(phase);
%                     Gesp = rho .* (1+phase) + conj(rho) .* phase - ifftnshift(kSpaceCombi,obj.trafo.fftdim,obj.trafo.scrambledim) .* (1+phase);
                end
            elseif(obj.espresso.constraint == 2) % min Im
                g = rho - abs(rho);
                if(obj.espresso.norm == 1) % L1
                    A = 0.5 .* conj(g)./(abs(g) + eps);
                    B = 1 - 0.5 .* conj(rho)./(abs(rho) + eps);
                    C = conj(A);
                    D = -0.5 .* conj(rho)./(abs(rho) + eps);
                    Gesp = conj(A.*B + C.*D);
                else % L2
                    f = sqrt(g(:)'*g(:));
                    A = conj(g)./(2*f + eps);
                    B = 1 - 0.5 .* conj(rho)./(abs(rho) + eps);
                    C = conj(A);
                    D = -0.5 .* conj(rho)./(abs(rho) + eps);
                    Gesp = conj(A.*B + C.*D);
                end
            end
        else
            Gesp = 0;
        end
        
        res = L2Grad + obj.lambda .* Gsparse + obj.lambdaTV .* Gtv + obj.lambdaESPReSSo .* Gesp;    
    end

    function res = objective(rho,dx,t)
        % L2 norm part
        w = kSpace - obj.fullMask.*fftnshift(rho+t*dx,obj.trafo.fftdim,obj.trafo.scrambledim);
        L2Obj = w(:)'*w(:);
        
        % L1 norm part
        % sparsify
        w = kernelFTrafo(rho+t*dx);
        L1Obj = sum((conj(w(:)).*w(:)+obj.l1Smooth).^(1/2));
        
        % TV
        if(obj.lambdaTV)
            w = rho+t*dx;
            Dx = w([2:end,end],:,:,:) - w;
            Dy = w(:,[2:end,end],:,:) - w;
            Dz = w(:,:,[2:end,end],:) - w;

            w = cat(5,Dx,Dy,Dz);
            TVObj = sum((conj(w(:)).*w(:)+obj.l1Smooth).^(1/2));
        else
            TVObj = 0;
        end
        
        % ESPReSSo
        if(obj.lambdaESPReSSo)
            if(obj.espresso.constraint == 1) % POCS
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
                w = rho+t*dx;
                rhoTmp = bTrafo(w);
                if(any(obj.trafo.fftBA(1,:))) % go back to kSpace
                    for iBA=find(obj.trafo.fftBA(1,:))
                        rhoTmp = fftnshift(rhoTmp,iBA);
                    end
                end
                [lmaskSym,pfDim,phase] = espresso_reconGrad(permute(rhoTmp,[4 1 3 2]), cs); % in: cha-k_y-k_x-k_z, out: y-z-x-cha (cart)
                phase = permute(phase,[1 3 2 4]);
                phase = exp(2i * angle(phase));
                lmaskSym = repmat(squeeze(lmaskSym),[1 1 1 nCha]);
%                 lmaskSym = permute(lmaskSym,[1 3 2 4]);
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

                rhoTmp = bTrafo(conj(w).*phase + w);
                if(any(obj.trafo.fftBA(1,:))) % go back to kSpace
                    for iBA=find(obj.trafo.fftBA(1,:))
                        rhoTmp = fftnshift(rhoTmp,iBA);
                    end
                end
                ESPObj = rhoTmp - kSpaceCombi;
                
            elseif(obj.espresso.constraint == 2) % Im(rho)
                w = rho+t*dx; 
                ESPObj = w - abs(w);
            end
                
            if(obj.espresso.norm == 1) % L1
                ESPObj = sum((conj(ESPObj(:)).*ESPObj(:)+obj.l1Smooth).^(1/2));
            else % L2
                ESPObj = ESPObj(:)'*ESPObj(:);
            end
        else
            ESPObj = 0;
        end
        
        
        res = L2Obj + obj.lambda .* L1Obj + obj.lambdaTV .* TVObj + obj.lambdaESPReSSo .* ESPObj;
    end
    
end

