function [smplPtrnSym, pfDim, phaseOut] = espresso_reconGrad( ksp, cs )
%ESPReSSo Reconstruction with CG -> calculate needed variables
%
% [im, kspFull] = pocs( kspIn, iter, watchProgr )
%
% === Input ===
%
%   kspIn:      Reduced Cartesian MRI Data-Set
%               Any dimension may be reduced,
%               but only one reduction dim. is allowed due to Physics/Math.
%
%               Allowed shapes for kspIn are...
%                 ... Ny x Nx
%                 ... Nc x Ny x Nx
%                 ... Nc x Ny x Nx x Nz
%
%               With Nc == number of receive Channels / Coils.
%
%               kspIn can either be a zero-padded array, so the partial Fourier property is obvious.
%               Or kspIn can be the measured data only, then we try to find k-space centre automagically
%               and create a zero-padded array with the full size, first.
%               Errors are however more likely to occur in the latter case.
%
%
%   iter:       No. of iterations
%   (optional)  default: iter = 20
%               Try on your own if larger iter improves your results!
%
%   watchProgr: true/false; Whether the progress of the reconstruction should
%   (optional)  be monitored in an image window.
%               In 3D data, only the central partition will be shown.
%
%
% === Output ===
%
%   im:         Reconstructed Images (channels not combined)
%
%   kspFull:    Reconstructed full k-space data
%
%
%
% === About the code ===
%
%   (1) We find out whether input data is
%       a) already zero-filled or
%       b) the pure asymmetric dataset, only
%
%       If b) is true, we zero-fill the data ourselves, which means we have to
%       determine the dimension first, in which the partial Fourier reduction was done.
%       We therefor find the position of the max. intensity in k-space which should
%       be identical to k-space centre. If the k-space centre is different from the
%       centre of the matrix, we know the partial Fourier dimension.
%       We then enlarge the matrix to its desired full size and fill the new part
%       with zeros.
%       If a) was true, finding the partial Fourier dimension is easy:
%       It is the dimension with all the zeros. :-)
%
%   (2) We create one low resolution image per channel/coil:
%
%       We need a symmetrically sampled part around the central k-space. Think of a
%       small stripe of phase encoding lines in the central k-space.
%       We only use these symmetric data (setting the rest zero) to reconstruct
%       low-resolution images. In order to avoid Gibbs-Ringing, a Hamming-filter
%       with the width of the stripe is multiplied with the data.
%       Additionally, all the fully sampled dimensions get a Hamming filter, too,
%       since we increase SNR, reduce further Gibbs-ringing and do not lose much
%       resolution.
%
%   (3) The phase of the low-resolution images is saved
%
%       POCS uses the fact that k-space data of real objects (no imaginary part)
%       have a point symmetry:
%           S(-k)  =  S*(k)      with k = (kx, ky, kz)
%       Our MRI objects are always complex, but we assume that phase variations
%       are due to coil sensitivities and B0-inhomgeneities,
%       which are both slowly varying (no high res. required).
%       Small-scale phase pertubations will decrease the reconstruction quality.
%
%   (4) Reference phase is applied in image space
%
%       We...
%       ... transform our zero-filled data to image space (IFFT)
%       ... remove the phase --> abs(image)
%       ... set the phase of our reference phase map --> image .* exp(1i.*phase)
%       ... transform back to k-space (FFT)
%       ... re-insert the measured data (self-consistency!)
%       ... goto "We..."
%
%       Iterating through the above steps fills the missing k-space points
%       with reasonable values.
%       If the phase varies slowly and there is no aliasing, this works very well.
%
% Aliasing artifacts are very challenging for POCS.
% So try to prevent aliasing in the first place (sufficient Field of View).

% =========================================================================
% Original code by Martin Blaimer
% * changed by Uvo Hoelscher
% * changed by Michael VÃ¶lker
%        -- auto-detect PF dimension
%        -- auto-find centre point/line/partition
%        -- accept zerofilled or "pure" data
%        -- for multichannel or plain 2D data (single-channel)
%        -- 2D and 3D
%        -- error handling
%        -- comments, comments, comments
%        -- added option to monitor progress
%        -- moved code to seperate functions
%        -- smooth transition between acquired signal and
%           reconstructed data
% * changed by Thomas Kuestner, University of Tuebingen
%   (thomas.kuestner@med.uni-tuebingen.de)
%        -- sparse sampling matrix for use in Compressed Sensing MRI -> espresso_recon
%

% =========================================================================

    % ( ===================================================================
    % Input Handling
    %
        if ~exist( 'ksp', 'var' ) || isempty(ksp) || ~isnumeric(ksp)
            error('pocs:input', 'First input must be Cartesian k-space data.')
        end
%         if ~exist('iter','var') || isempty(iter) || numel(iter) ~= 1 || ~isnumeric(iter)
%             iter = 20;
%         end
%         if ~exist('watchProgress','var')  || isempty(watchProgress) ||  numel(watchProgress) ~= 1 || ~isfinite(watchProgress)
%             watchProgress = false;
%         else
%             watchProgress = logical( watchProgress );
%         end
        
        % TK:
        if( ~exist('cs','var') || isempty(cs) )
            nupf = false; % non-uniform partial fourier
%             cs.skipDetection = false;
        else
            nupf = true;
        end

        Ndim = ndims( ksp );

        if Ndim > 4 || Ndim < 2
            error('pocs:shape','First input ''kspace'' should have one of these shapes:\n\n\t... Ny x Nx\n\t... Nc x Ny x Nx\n\t... Nc x Ny x Nx x Nz')
        end
        if Ndim == 2    % Ny x Nx
            ksp = reshape( ksp, [1 size(ksp)] );    %  1 x Ny x Nx  --> now we have one channel...
            wasAddedCoilDim = true;
            Ndim = 3;
        else
            wasAddedCoilDim = false;
        end

        % read the properties of the data
        sz   = size( ksp );
        sz   = sz(2:end);           % the (k-)spatial size of the array (i.e. without channels)
        prec = class( ksp );        % single or double precision?
    % ) ===================================================================

    % First: Check the sampling pattern (which parts of input are actually data?)
    smplPtrn = reshape( sum(abs(ksp),1) ~= 0, sz);        % Ny x Nx x Nz


    % ( ===================================================================
    % If input data is not yet zero-filled, do it here
    %
    % TK:
    if(~nupf)
        if nnz(smplPtrn) == numel(smplPtrn)     % only the sampled data were passed / |N|umber of |N|on |Z|ero elements

            [ ksp, pfDim, isUpper, isLower, Nsmp ] = zerofillPFdim( ksp, wasAddedCoilDim );

            sz = size( ksp );
            sz = sz(2:end);     % ignore channels
        else
            [ pfDim, isUpper, isLower, Nsmp ] = detectPFdim( smplPtrn, wasAddedCoilDim );
        end
    else
%         pfDim = cs.pfDim;
%         isUpper = cs.isUpper;
%         isLower = cs.isLower;
%         Nsmp = cs.Nsmp;
        % smplPtrn: y-x-z or y-x
        smplPtrn = cs.smplPtrn; % sparse sampling pattern before 1st CS recon
        pfn = cs.pfn;
        % out: smplPtrn: smooth mask (complete 0/1 regions)
        [ smplPtrn, smplPtrnSparse, pfDim, isUpper, isLower ] = detectPFdimCS(smplPtrn, wasAddedCoilDim, pfn);
%         smplPtrnSparse = cs.smplPtrnSparse; % sparse sampling pattern of current CS recon loop
    end

% TK: %         clear  smplPtrn
    % ) ===================================================================


    if numel(sz) < 3
        sz(3) = 1;
    end
    Ny = sz(1);
    Nx = sz(2);
    Nz = sz(3);


    % ( ===================================================================
    % Handle ugly problems.
    %
        if ~isUpper && ~isLower
            error('pocs:UnknownErrorFound', 'I thought we are partial Fourier, but things seem to make no sense... :-(')
        end
    % ) ===================================================================



    % =====================================================================
    %
    %        We can now be sure to operate with zero-padded data.
    %
    % =====================================================================



    % initialize a cell of subscripts
    subs = { ':', ':', ':', ':' };      % all channels / all Ny / all Nx / all Nz

    % If the first entries are zero-filled (instead of the trailing ones),
    % flip the entries so we can treat them as if we pf'ed the first half of kspace.
    if isLower
        subs{pfDim+1} = sz(pfDim):-1:1;     % ...esreveR
        ksp           = ksp(subs{:});       % !ecaps-k si sihT
        smplPtrn      = smplPtrn(subs{2:end});
        if(nupf)
            smplPtrnSparse= smplPtrnSparse(subs{2:end});
        end
%         if(nupf)
%             smplPtrn  = smplPtrn(subs{2:end});
%         else
%             smplPtrn  = smplPtrn(subs{2:end});  % TK
%         end
        subs{pfDim+1} = 1:sz(pfDim);        % lalala, we didn't do anything...
    end

    % TK
    if(nupf)
        [applDim,idx,~] = getCalibrationSize(ksp, smplPtrn, pfDim); % Ny x Nx x Nz
    else
        % Find out which point is in the centre and which indices belong to the
        % symmetrically sampled part of k-space.
        [ centreLine, idxSym ] = findSymSampled( ksp, pfDim, Nsmp );
        
        szSym = numel( idxSym );                % 2 * (Nsmp - centreLine) + 1

%         if isUpper
%             fprintf('Using %g points around point %g\n', szSym,    centreLine    );
%         else
%             fprintf('Using %g points around point %g\n', szSym, sz(pfDim)-centreLine+1 );
%         end
    
    end

    % ( ===================================================================
    % build up a symmetric low-pass filter
    %
    % TK
    if(nupf)
        filter = cast(ones(size(ksp)),prec);
    else
        filter = cast( 1, prec );
    end
    for d = 1:Ndim-1
        
        % TK
        if(nupf)
            reshRule = ones(1,Ndim);    % how the filter will be reshaped
            repRule = size(ksp);
            
            if d ~= pfDim   % Each standard dimension gets a simple low-pass filter

                filt1D = hamming( sz(d), 'periodic'  );
                
                % reshape the filter according to the dimension it represents
                reshRule(d+1) = sz(d);
                repRule(d+1) = 1;
                
                filt1D = reshape(filt1D, reshRule);
     
            else
                % create a narrow filter and remove everything else
                filt1D = zeros(sz(d), sz(applDim == 2), prec);            % full-size filter
                for l=1:sz(applDim == 2)
                    tmp              = hann( numel(idx{l}) + 2, 'symmetric' ); 
                    filt1D(idx{l},l) = tmp(2:end-1);
                end
                
                % reshape the filter according to the dimension it represents
                loopDim = find(applDim == 2);
                repRule(d+1) = 1;
                repRule(loopDim + 1) = 1;
                permRule = zeros(1,Ndim);
                permRule(find(applDim == 1)+1) = 1;
                permRule(find(applDim == 2)+1) = 2;
                permRule(permRule == 0) = 3:Ndim;
                
                filt1D = permute(filt1D,permRule);
%                 filt1D = shiftdim(filt1D,-d);
            end
                        
            filt1D = repmat(filt1D, repRule);
                  
            filter = filter .* filt1D;
                
        else
            reshRule = ones(1,Ndim);    % how the filter will be reshaped

            if d ~= pfDim   % Each standard dimension gets a simple low-pass filter

                filt1D = hamming( sz(d), 'periodic'  );

            else            % our partial Fourier dimension gets an extra nice filter

                % create a narrow filter and remove everything else
                filt1D          = zeros(sz(d), 1, prec);            % full-size filter
                tmp             = hann( szSym + 2, 'symmetric' );   % a very narrow window
                filt1D(idxSym)  = tmp(2:end-1);                     % cut out the zeros at the edges (we have data there!)
               
                % take a look:
                %figure, plot(filt1D)
            end

            % reshape the filter according to the dimension it represents
            reshRule(d+1) = sz(d);

            filt1D = reshape( filt1D, reshRule );
            filter = bsxfun( @times, filter, filt1D );      % iteratively build up a multidimensional filter
        end
    end
    % ) ===================================================================

    % Apply the low-pass filter
    if(nupf)
        kspLowRes = ksp .* filter;
    else
        kspLowRes = bsxfun( @times, filter, ksp);
    end
    clear  filt1D  filter  reshRule  idxSym


    % ( ===================================================================
    % prerequisites prior to the iteration loop
    %
    %  Set everything up here, do computations that you don't have
    %  to do in the loop, remove no longer needed variables...
    %
        % fftshift everything once before and after for-looping
        %  => less overhead during iteration
        %ksp       = cmshiftnd(       ksp, [0  sz/2] );
        kspLowRes = cmshiftnd( kspLowRes, [0  sz/2] );

        % reorder arrays such that the fft-dimensions come first
        % => faster memory access
        %ksp       = permute(       ksp, [2 3 4 1] );      % Ny x Nx x Nz x Nc
        kspLowRes = permute( kspLowRes, [2 3 4 1] );      %
        subs      = { subs{2}, subs{3}, subs{4}, subs{1} };

        % calc. initial image and the reference phase map
        %im        =  fft( fft( fft( conj(ksp), [], 1), [], 2), [], 3);  % im's phase is wrong now, but we only want it's abs() to be correct
        phaseOut  = ifft(ifft(ifft( kspLowRes, [], 1), [], 2), [], 3);
        %phase     = exp(1i * angle(phaseOut));
%         phase     = phaseOut./abs(phaseOut);

        % We use a trick in the loop to avoid using ifft (fft is faster).
        % We only need to calculate the factor 1/N ourselves, with N = prod(sz)
        %phase = phase ./ prod(sz);      % 1/N is absorbed inside the phase array, once
        
        % create image with calculated phasemap from low res image
        %im = abs(im) .* phase;

        % In the loop, we want to know where we have to copy the
        % measured data to, so we set the subscript of the pf dimension
        % accordingly.
        % We have to do this due to the ifftshift'ing above.
%         if(nupf)
%             smplPtrnSparse = repmat(smplPtrnSparse, [ones(1,length(sz)), size(ksp,4)]);
%             smplPtrnSparse = cmshiftnd(smplPtrnSparse, [sz/2 0]);  
%             smplPtrn = repmat(smplPtrn, [ones(1,length(sz)), size(ksp,4)]);
%             smplPtrn = cmshiftnd(smplPtrn, [sz/2 0]);            
%             tmp         = true( 1, sz(pfDim));
%             subs{pfDim} = find(ifftshift(tmp));
%         else
%             tmp         = false( 1, sz(pfDim));
%             tmp(1:Nsmp) = true;
%             subs{pfDim} = find(ifftshift(tmp));
%         end
%         
%         % release RAM
%         clear  tmp  kspLowRes
% 
%         % only keep the acquired data in memory
%         if(~nupf)
%             ksp = ksp(subs{:});
%         end
%     % ) ===================================================================

    % Helpers for pretty-printing:
    % Such a mess for such beautiful output!
%     dispProgress('ESPReSSo', 0, iter); 
%     b = repmat('=',1,80);
%     progress_str = 'starting POCS loop...';
%     fprintf( '%s\n%s\n%s   %s', b, b(1), b(1), progress_str )
%     edging = sprintf( '\n%s\n%s', b(1), b );
%     fprintf( edging )


    % ( ===================================================================
    % iterative reconstruction POCS
    %
%     tic
%     for ii = double(~watchProgress) : iter
% 
%         if ii > 0
% 
%             % Fourier transform the image to k-space
%             im = fft(fft2(im),[],3);
% %             im = fft(fft(fft(  im  ,[],1),[],2),[],3);      % "im" is a really bad variable name now
%                                                             % but we save a lot of RAM with this
%             % Data Consistency:
%             % insert original data where we have them
%             if(nupf)
%                 im(smplPtrn) = ksp(smplPtrn);
%             else
%                 im(subs{:}) = ksp;                              % "im" is still our reconstructed k-space signal
%             end
% 
%             % Fourier transform into image domain
%             im = conj( im );
%             im = fft(fft2(im),[],3);
% %             im = fft(fft(fft(  im  ,[],1),[],2),[],3);      % Now, "im" is an image again.
% 
%             % create image with calculated phasemap from low res image
%             im = abs(im) .* phase;
%     
%             dispProgress('ESPReSSo', ii/iter);
% %             prevLength = numel(progress_str) + numel(edging);
% %             t = toc;
% %             ETA = (t./ii) * iter  - t;
% %             progress_str = sprintf( 'Iteration %g/%g, in %g s,  ETA: %g s...', ii, iter, t, ETA );
% % 
% %             fprintf([repmat('\b',1,prevLength) '%s' '%s'], progress_str, edging );
% 
%         end % if ii > 0
% 
%         % a rough way to monitor the progress
%         %
%         if watchProgress
%             tmp = ifftshift(sqrt(sum(abs(im(:,:,1,:).^2),4)));      % due to fftshift(), the 1st partition is the central one
%             maxRange = sort( tmp(:), 'descend' );
%             maxRange = maxRange( ceil(0.05 * numel(maxRange)) );    % ignore the "hottest" 5%
% 
%             if ~exist('pic','var')
%                 pic = [tmp tmp zeros(size(tmp),prec)];
%                 diffScale = 1;
%             else
%                 delta = abs( pic(:,Nx+(1:Nx)) - tmp );
%                 diffScale = 0.5 * maxRange / median( delta(:) );
%                 pic(:,  Nx+(1:Nx)) = tmp;
%                 pic(:,2*Nx+(1:Nx)) = diffScale * delta;
%                 clear delta
%             end
% 
%             figure(999)
%             imagesc( pic, [0    maxRange ] )
%             title(sprintf('\\bfiteration %g\ninitial     |    current     |     abs(previous - current) — %g', ii, diffScale ))
%             axis image
%             colormap(gray(256))
%             drawnow
%             clear tmp
% 
%             %if Nz == 1          % little pause for 2D (too fast otherwise)
%             %    pause(2 / iter)
%             %end
%         end
% 
%     end % for ii = 1:iter
% %     fprintf([repmat('\b',1, numel(progress_str) + numel(edging)) 'POCS done! (%g s)' '%s\n\n'], t, edging );
%     dispProgress('ESPReSSo', 'Close');
    % ) ===================================================================

%     clear pic % phase
    

    % ( ===================================================================
    % The main part is over. Time for some thoughts.
    %
    % We began with a dataset that had fewer data samples than would be
    % necessary for an unambiguous image reconstruction. As a consequence,
    % an infinite number of images corresponds to the acquired data.
    % The above iteration picks that single image whose abs() fits the data
    % AND whose phase corresponds to the low-resolution phase, obtained
    % using the symmetric part of the data.
    %
    % Viewed in k-space, there is almost always a severe edge at the border
    % between acquired and interpolated data, which is due to imperfections
    % in the assumptions made.
    % Namely, phase often has some high frequency components which cannot be
    % accounted for in the low-resolution map. Additionally, there is noise
    % and we may have changing contrast or trajectory errors in our MRI
    % sequence.
    %
    %         ^
    %         | A A A A A A A A   \
    %         | A A A A A A A A
    %         | A A A A A A A A     acquired signal
    %      k2 | A A A A A A A A
    %         | A A A A A A A A   /
    %         | I I I I I I I I   \
    %         | I I I I I I I I     interpolated data
    %         | I I I I I I I I   /
    %         ----------------->
    %                 k1
    %
    % Empirically, it should be wise to create a smoother transition from
    % the acquired part of the signal to the interpolated data.
    %
        if(nupf)
%             loopDim = find(applDim == 2);
% %             im = fft(fft(fft(  im  ,[],1),[],2),[],3);      % "im" becomes k-space signal, again
%             for l=1:sz(applDim == 2)
% %                 Ntrans = floor((numel(idx{l})-1)/2);
%                 Ntrans = numel(idx{l});
%                 Nsmp   = idxBorder(1,l);
%                 
%                 % Create subscripts where we intend to keep the measured data, only.
%                 tmp                 = false( 1, sz(pfDim));
%                 tmp(1:Nsmp-Ntrans)  = true;
%                 subsPure            = subs;
%                 subsPure{pfDim}     = find(ifftshift(tmp));
%                 subsPure{loopDim}   = l;
%                 
%                 % Create subscripts where we want to have a smooth transition between
%                 % measured and phase-corrected data.
%                 subsTrans           = subs;
%                 tmp                 = false( 1, sz(pfDim));
%                 tmp(idx{l})         = true;
% %                 subsTrans{pfDim}    = setdiff( subs{pfDim}, subsPure{pfDim} );
%                 subsTrans{pfDim}    = find(ifftshift(tmp));
%                 subsTrans{loopDim}  = l;
%                 
%                 % build a filter for the transition:
%                 tmp         = hann( Ntrans+2, 'symmetric');
%                 filterTrans = tmp( 2 : end-1 );
%                 filterTrans = reshape( filterTrans, [ ones(1,pfDim-1)  Ntrans  1] );
%                 
%                 % Seperate data in unfiltered part and transition zone.
% %                 tmp = zeros( size(im), prec );
% %                 tmp(subs{:}) = ksp;             % kspace domain
%                 kspPure  = ksp(subsPure{:});
%                 kspTrans = ksp(subsTrans{:});
%                 
%                 
%                 im(subsPure{:}) = kspPure;                      % strict data consistency for Nsmp-Ntrans samples
%                 im(subsTrans{:}) =  bsxfun( @times,   filterTrans,  kspTrans         )     ...
%                                   + bsxfun( @times, 1-filterTrans,  im(subsTrans{:}) );   
%             end
%             
%             if nargout > 1
%                 kspFull = im;
%             else
%                 kspFull = double.empty([sz 0]);     % kspFull exists, but no memory required
%             end
%             im = ifft(ifft(ifft(  im  ,[],1),[],2),[],3);   % "im" is an image, again
            
            % put out a symmetrically sampled mask
            smplPtrnSym = smplPtrnSparse;
%             if(pfDim == 1)
%                 smplPtrnSym = smplPtrnSparse | smplPtrnSparse(end:-1:1,:,:,:);
%             elseif(pfDim == 2)
%                 smplPtrnSym = smplPtrnSparse | smplPtrnSparse(:,end:-1:1,:,:);
%             elseif(pfDim == 3)
%                 smplPtrnSym = smplPtrnSparse | smplPtrnSparse(:,:,end:-1:1,:);
%             end
            
        else
%             Ntrans = floor( (szSym-1)/3 );      % width of the transition zone
% 
%             % Create subscripts where we intend to keep the measured data, only.
%             tmp                 = false( 1, sz(pfDim));
%             tmp(1:Nsmp-Ntrans)  = true;
%             subsPure            = subs;
%             subsPure{pfDim}     = find(ifftshift(tmp));
% 
%             % Create subscripts where we want to have a smooth transition between
%             % measured and phase-corrected data.
%             subsTrans           = subs;
%             subsTrans{pfDim}    = setdiff( subs{pfDim}, subsPure{pfDim} );
% 
%             % build a filter for the transition:
%             tmp         = hann( 2*Ntrans+3, 'symmetric');
%             filterTrans = tmp( Ntrans+3 : end-1 );
%             filterTrans = reshape( filterTrans, [ ones(1,pfDim-1)  Ntrans  1] );
% 
%             % Seperate data in unfiltered part and transition zone.
%             tmp = zeros( size(im), prec );
%             tmp(subs{:}) = ksp;             % kspace domain
%             kspPure  = tmp(subsPure{:});
%             kspTrans = tmp(subsTrans{:});
%             clear  tmp  ksp
% 
%             im = fft(fft(fft(  im  ,[],1),[],2),[],3);      % "im" becomes k-space signal, again
% 
%             im(subsPure{:}) = kspPure;                      % strict data consistency for Nsmp-Ntrans samples
%             im(subsTrans{:}) =  bsxfun( @times,   filterTrans,  kspTrans         )     ...
%                               + bsxfun( @times, 1-filterTrans,  im(subsTrans{:}) );
% 
%             clear  subsPure  subsTrans  filterTrans  kspPure  kspTrans
% 
%             if nargout > 1
%                 kspFull = im;
%             else
%                 kspFull = double.empty([sz 0]);     % kspFull exists, but no memory required
%             end
% 
%             im = ifft(ifft(ifft(  im  ,[],1),[],2),[],3);   % "im" is an image, again
            
            % put out a symmetrically sampled mask
            smplPtrnSym = [];
        end
    % ) ===================================================================

    

    % ( ===================================================================
    % Undo the prerequisites (--> postrequisites???)
    %
        % undo the permutations
%         im      = permute(      im, [4 1 2 3] );
%         kspFull = permute( kspFull, [4 1 2 3] );
        subs    = { subs{4}, subs{1}, subs{2}, subs{3} };

        % undo the fftshifts
%         im      = cmshiftnd(      im, [0  sz/2] );
%         kspFull = cmshiftnd( kspFull, [0  sz/2] );
        
        if(~isempty(smplPtrnSym))
            smplPtrnSym = permute( smplPtrnSym, [4 1 2 3] );
%             smplPtrnSym = cmshiftnd( smplPtrnSym, [ 0 sz/2] );
        end
        % undo flipping
        if isLower
            subs{pfDim+1} = sz(pfDim):-1:1;
%             im            = im(subs{:});
%             kspFull       = kspFull(subs{:});
            if(~isempty(smplPtrnSym))
                smplPtrnSym = smplPtrnSym(subs{:});
                if(pfDim == 1)
                    smplPtrnSym = smplPtrnSym | smplPtrnSym(:,end:-1:1,:,:);
                elseif(pfDim == 2)
                    smplPtrnSym = smplPtrnSym | smplPtrnSym(:,:,end:-1:1,:);
                elseif(pfDim == 3)
                    smplPtrnSym = smplPtrnSym | smplPtrnSym(:,:,:,end:-1:1);
                end
            end
        end
    % ) ===================================================================

%     if wasAddedCoilDim                              % we initially reshaped a simple 2D raw data matrix to be of size 1 x Ny x Nx
%         im      = reshape(      im, Ny, Nx, [] );
%         kspFull = reshape( kspFull, Ny, Nx, [] );
%     end

end     % of pocs()






% =========================================================================
%                                                                         =
%                      SWAPPED CODE                                       =
%                                                                         =
% =========================================================================






function [ ksp, pfDim, isUpper, isLower, Nsmp ] = zerofillPFdim( ksp, wasAddedCoilDim )
    % Only the acquired data were passed and we have to find the asymmetric
    % dimension. Then we increase the size along this dimension and pad with 0.

    Ndim    = ndims( ksp ) - 1;     % one dimension was for the channels
    sz      = size(  ksp );
    sz      = sz(  2:end );         % ignore channel dimension
    Nc      = size(  ksp, 1 );
    prec    = class( ksp );

    % init some helper variables
    pfDim    = 0;               % partial Fourier reduction dimension
    isUpper  = false;
    isLower  = false;
    isPartialFourier = false(Ndim,1);

    % ( ===============================================================
    % autodetect the Partial Fourier dimension
    %
    for d = 1:Ndim

        centre = floor( sz(d)/2 ) + 1;

        tmp = squeeze( sum(abs(ksp),1) );
        for d2 = 1:Ndim
            if d2 ~= d
                tmp = max(tmp,[],d2);       % keep only the maximum of non-partial data points
            end
        end
        [ dummy, maxPos(d) ] = max( tmp(:) );   %#ok <-- don't use "~", for compatibility

        if abs(maxPos(d) - centre) >= 2      % significant asymmetry ==> partial Fourier acquisition

            isPartialFourier(d) = true;
            pfDim = d;
            Nsmp = sz(d);

            isUpper = maxPos(d) > centre;   % Did we sample the upper matrix part, so the lower part is missing...
            isLower = maxPos(d) < centre;   % ... or are the first data points missing (e.g. asymmetric echo)?
        end
    end % for d = 1:Ndim
    %
    % ) ===== (PF dim detection) ======================================


    switch nnz(isPartialFourier)    % |N|umber of |N|on |Z|ero elements
        case 0
            error( 'pocs:NoPfDim', 'No partial Fourier dimension found.' )
        case 1
%             fprintf( 'Found partial Fourier along array dimension %d\n', pfDim + ~wasAddedCoilDim )
        otherwise
            error( 'pocs:TooManyPfDims', 'Partial Fourier only allowed in 1 dimension, but %g were found!', nnz(isPartialFourier) )
    end

    if pfDim == 0   % our init value above
        error('zerofillPF:NoPF','No partial Fourier property found!')
    end

    % initialize a cell of subscripts
    subs = { ':', ':', ':', ':' };      % all channels / all Ny / all Nx / all Nz

    c = maxPos(pfDim);

    if isUpper

        sz(pfDim) = 2 * (c - mod(c,2));         % determine the blown-up size we want to achieve
        subs{pfDim+1} = 1:Nsmp;

    elseif  isLower

        sz(pfDim) = 2 * (Nsmp - c + 1);
        c = floor( sz(pfDim)/2 ) + 1;
        sz(pfDim) = sz(pfDim) + 2*~mod(c,2);    % A hack for Stefan's data... keep an eye on this!
        subs{pfDim+1} = (1:Nsmp) + (sz(pfDim)-Nsmp);

    else
        error( 'zerofillPF:PFdimNotClassified', 'Could not tell how partial Fourier was implemented.' )
    end

    % do the zerofilling
    tmp = zeros( [Nc sz], prec );
    tmp(subs{:}) = ksp;
    ksp = tmp;

end % of zerofillPFdim()

function [ pfDim, isUpper, isLower, Nsmp ] = detectPFdim( smplPtrn, wasAddedCoilDim )
    % User passed already zero-padded data. This was nice, now it's easy
    % to find the partial Fourier dimension!

    Ndim = ndims( smplPtrn );
    sz = size( smplPtrn );

    % init some helper variables
    pfDim    = 0;               % partial Fourier reduction dimension
    isUpper  = false;
    isLower  = false;
    isPartialFourier = false( Ndim, 1 );

    % ( ===============================================================
    % Determine if this is a zerofilled partial Fourier measurement
    % and along which dimension the data is reduced.
    %
    % smplPtrn in Partial Fourier looks like this:
    %
    %     ^
    %     | 1 1 1 1 1 1 1 1     --->  sampling pattern is the same
    %     | 1 1 1 1 1 1 1 1           for all k1 points
    %     | 1 1 1 1 1 1 1 1
    %  k2 | 1 1 1 1 1 1 1 1             i.e. for programming:
    %     | 1 1 1 1 1 1 1 1             smplPtrn == repmat( smplPtrn(:,1,1), [1 Nx Nz] )
    %     | O O O O O 0 0 0
    %     | O O O O O 0 0 0
    %     | O O O O O 0 0 0
    %      ----------------->
    %           k1
    %
    for d = 1:Ndim

        subs = { ones(1,sz(d)),     ... % initialize a cell of subscripts we might be interested in
                 ones(1,sz(d)),     ... % TK: subs = cell(1,Ndim); subs = deal(ones(1,sz(d));
                 ones(1,sz(d))  };
        subs{d} = 1:sz(d);              % we ask for all entries in the d'th dimension

        idx_d = sub2ind( sz, subs{:} ); % convert to linear array indices

        oneCol = smplPtrn( idx_d );     % one column of the d'th dimension

        % create a rule how to reshape oneCol
        reshRule = ones(1,Ndim);
        reshRule(d) = sz(d);                    % e.g. reshRule = [   1   1 128 ]
        oneCol = reshape( oneCol, reshRule);

        % create a rule how to replicate oneCol
        repRule = sz;
        repRule(d) = 1;                         % e.g. repRule  = [ 256 256   1 ]

        % Check if we get the sampling pattern again
        % just by replicating oneCol along the other dimensions
        isPartialFourier(d) = isequal( smplPtrn, repmat( oneCol, repRule ) );

        if isPartialFourier(d)
            pfDim = d;
            Nsmp  = nnz( oneCol );      % how many fully sampled lines do we have?

            % Sampled upper or lower part of k-space matrix?
            isUpper = isequal( oneCol(:).', [ true( 1,Nsmp)          false(1,sz(d)-Nsmp)    ]);
            isLower = isequal( oneCol(:).', [ false(1,sz(d)-Nsmp)    true( 1,Nsmp)          ]);
        end
    end
    % ) ===============================================================

    switch nnz(isPartialFourier)    % |N|umber of |N|on |Z|ero elements
        case 0
            error( 'pocs:NoPfDim', 'No partial Fourier dimension found.' )
        case 1
%             fprintf( 'Found partial Fourier along array dimension %d\n', pfDim + ~wasAddedCoilDim )
        otherwise
            error( 'pocs:TooManyPfDims', 'Partial Fourier only allowed in 1 dimension!' )
    end

end % of detectPFdim()

function [ centreLine, idxSym ] = findSymSampled( ksp, pfDim, Nsmp )

    Ndim = ndims( ksp ) - 1;    % one for channels
    sz = size( ksp );
    sz = sz(2:end);

    % autodetect the central k-space line
    %if ~exist('centreLine', 'var') || isempty(centreLine)
        tmp = squeeze( sum(abs(ksp),1) );
        for d = 1:Ndim
           if d ~= pfDim
               tmp = max(tmp,[],d);     % keep only the maximum of non-partial data points
           end
        end
        [ dummy, centreLine] = max( tmp(:) );   %#ok the central line has the max intensity
    %end

    % calculate the size of the symmetric part and the full dataset
    startSym = centreLine - (Nsmp - centreLine);    % start of our symmetric sampling
    endSym   = centreLine + (Nsmp - centreLine);    % end of symmetric part
    idxSym   = startSym : endSym;

    if any(idxSym < 1)    ||   any(idxSym > sz(pfDim))
       error( 'pocs:BadDataProperty' , 'Symmetric part of k-space out of bounds.\nThe maximum k-space intensity is at index %g whereas it should be centred => near %g.\nThe way, zerofilling was done is probably wrong.\nCheck your input k-space.', centreLine, round(sz(pfDim)/2) )
    end

end % of findSymmetricSampled()

function x = cmshiftnd( x, shifts)
%Function to circularly shift N-D arrays

    if nargin < 2 || all(shifts(:) == 0)
       return                       % no shift
    end

    sz      = size( x );
    numDims = ndims(x);             % number of dimensions
    idx = cell(1, numDims);         % creates cell array of empty matrices,
                                    % one cell for each dimension

    for k = 1:numDims

        m = sz(k);
        p = ceil(shifts(k));

        if p < 0
            p = m + p;
        end

        idx{k} = [p+1:m  1:p];
    end

    % Use comma-separated list syntax for N-D indexing.
    x = x(idx{:});

end % of cmshiftnd()

% Avoid the need for the signal toolbox and implement
% hamming() and hann() manually:
%
function w = hamming( N, symFlag )
%Hamming window
%
% w = hamming(L) returns an L-point symmetric Hamming window in the column vector w.
% L should be a positive integer.
%
%  The coefficients of a Hamming window are computed from the following equation:
%
%       w(n) = 0.54  +  0.46 * cos(2*pi*n/N),   0 <= n <= N
%
%
% w = hamming( L, 'symFlag') returns an L-point Hamming window using the window sampling
% specified by 'symFlag', which can be either 'periodic'  or 'symmetric' (the default).
% The 'periodic' flag is useful for DFT/FFT purposes, such as in spectral analysis.
% The DFT/FFT contains an implicit periodic extension and the periodic flag enables a signal
% windowed with a periodic window to have perfect periodic extension.
% When 'periodic' is specified, hamming computes a length L+1 window and returns the first L points.
% When using windows for filter design, the 'symmetric' flag should be used.
%
% --> http://www.mathworks.de/de/help/signal/ref/hamming.html
% --> https://de.wikipedia.org/wiki/Hamming-Fenster

% implemented by Michael.Voelker@mr-bavaria.de, 2012

    if ~exist( 'N', 'var' ) || isempty(N) || numel(N) ~= 1 || ~isnumeric(N)  || ~isfinite(N) || N < 1 || floor(N) ~= N
        error( 'hamming:badSize', 'Window lenght must be a positive integer.' )
    end
    if ~exist( 'symFlag', 'var' ) || isempty(symFlag)
        symFlag = 'symmetric';
    end

    if N == 1
        w = 1;
        return
    end

    switch symFlag
        case 'symmetric'
            L = N-1;
        case 'periodic'
            L = N;
        otherwise
            error('hamming:symFlag', 'Unknown symmetry flag. Try ''symmetric'' (default) or ''periodic''.')
    end

    w = (0:N-1) - L/2;
    w = 0.54  +  0.46 * cos(2*pi * w(:)./L);

end % of hamming()


function w = hann( N, symFlag )
%von-Hann (Hanning) window
%
% w = hann(L) returns an L-point symmetric Hann window in the column vector w.
% L must be a positive integer.
%
% The coefficients of a Hann window are computed from the following equation:
%
%      w(n) = 0.5 * (1 + cos(2*pi*n/N)),   0 <= n <= N
%
% The window length is L = N+1.
%
% w = hann(L,'sflag') returns an L-point Hann window using the window sampling specified by 'sflag',
% which can be either 'periodic' or 'symmetric' (the default). The 'periodic' flag is useful for DFT/FFT purposes,
% such as in spectral analysis.
% The DFT/FFT contains an implicit periodic extension and the periodic flag enables a signal windowed
% with a periodic window to have perfect periodic extension.
% When 'periodic' is specified, hann computes a length L+1 window and returns the first L points.
% When using windows for filter design, the 'symmetric' flag should be used.
%
% --> http://www.mathworks.de/de/help/signal/ref/hann.html
% --> https://de.wikipedia.org/wiki/Hann-Fenster

% implemented by Michael.Voelker@mr-bavaria.de, 2012

    if ~exist( 'N', 'var' ) || isempty(N) || numel(N) ~= 1 || ~isnumeric(N)  || ~isfinite(N) || N < 1 || floor(N) ~= N
        error( 'hann:badSize', 'Window lenght must be a positive integer.' )
    end
    if ~exist( 'symFlag', 'var' ) || isempty(symFlag)
        symFlag = 'symmetric';
    end

    if N == 1
        w = 1;
        return
    end

    switch symFlag
        case 'symmetric'
            L = N-1;
        case 'periodic'
            L = N;
        otherwise
            error('hann:symFlag', 'Unknown symmetry flag. Try ''symmetric'' (default) or ''periodic''.')
    end

    w = (0:N-1) - L/2;
    w = 0.5 * ( 1 + cos(2*pi * w(:)./L) );

end % of hann()


% TK:

function [ smplPtrn, smplPtrnSparse, pfDim, isUpper, isLower ] = detectPFdimCS(mask, wasAddedCoilDim, pfn, type)
% find out the partial fourier dimension
% pfDim: 1 = y, 2 = x, 3 = z

Ndim = ndims(mask);
sz = size(mask); % y-x-z

if(nargin < 4)
    type = 'sparse';
end

if(strcmp(type,'sparse'))
    % sparse sampling pattern

    if(Ndim == 2)
        smplPtrn = mask;
        smplPtrnSparse = mask;

        % project onto rows (y-axis)
        projY = unique(sum(mask,2));
%         % project onto columns (x-axis)
%         projX = unique(sum(mask,1));
        if(length(projY) == 2)
            pfDim = 1;
            if(mod(size(mask,1),2) == 0)
                cut = [size(mask,1)/2, size(mask,1)/2+1];
            else
                cut = [ceil(size(mask,1)/2), ceil(size(mask,1)/2)];
            end
            halfside = [sum(mask(1:cut(1),:),2), sum(mask(cut(2):end,:),2)];
        else
            pfDim = 2;
            if(mod(size(mask,2),2) == 0)
                cut = [size(mask,2)/2, size(mask,2)/2+1];
            else
                cut = [ceil(size(mask,2)/2), ceil(size(mask,2)/2)];
            end
            halfside = [sum(mask(:,1:cut(1)),1), sum(mask(:,cut(2)+1:end),1)];
        end
    
        tmpDim = [2 1];
        if(halfside(1) > halfside(2))
            isUpper = true;
            isLower = false;
            idxBorder = find(sum(mask,tmpDim(pfDim)),1,'last'); % for y (pfDim=1) sum up over the rows (tmpDim=2)
            if(pfDim == 1)
                smplPtrn(1:idxBorder,:) = true;
            elseif(pfDim == 2)
                smplPtrn(:,1:idxBorder) = true;
            end
        else
            isUpper = false;
            isLower = true;
            idxBorder = find(sum(mask,tmpDim(pfDim)),1,'first');
            if(pfDim == 1)
                smplPtrn(idxBorder:end,:) = true;
            elseif(pfDim == 2)
                smplPtrn(:,idxBorder:end) = true;
            end
        end        

    else % Ndim >= 3
    
        sz = [sz, sz];
        dim = [1:Ndim, 1:Ndim];

        bordermap = cell(1,Ndim);

        for d=1:Ndim
            subs = {':',':',':'};

            % bordermap{d} is projection onto orthogonal dimensions of d containing
            % first nonzero element (border)
            % e.g.: d := y  => projection onto x-z plane and entries indicate first nonzero point in y direction 
            bordermap{d} = zeros(sz(d+1),sz(d+2));
            for e=1:size(bordermap{d},1)
                subs{dim(d+1)} = e;
                for f=1:size(bordermap{d},2)
                    subs{dim(d+2)} = f;
                    tmp = find(mask(subs{:}), 1, 'first');
                    if(~isempty(tmp))
                        bordermap{d}(e,f) = tmp;
                    end
                end         
            end   
        end

        % extract original 2D pattern
        % mask is the one which contains just "0" and "1" in the projected data
        helper = cellfun(@(x) unique(x(:)), bordermap, 'UniformOutput', false);
        helperVal = cellfun(@(x) ismember(0,x) & ismember(1,x), helper); 
        helperLen = cellfun(@(x) length(x), helper);
        tmpPos = helperLen == 2 & helperVal;
        mask2D = bordermap{tmpPos};
        idx_plane = find(tmpPos);

        if(idx_plane == 1)
            % plane 1: x-z plane
            pfDimIn = [2 3];
        elseif(idx_plane == 2)
            % plane 2: z-y plane
            pfDimIn = [3 1];
        elseif(idx_plane == 3)
            % plane 3: y-x plane
            pfDimIn = [1 2];
        else
            error('partialFourier(): Unknown projection plane');
        end

        % amount of sampled points (1st plane) and zero elements (2nd plane) per quadrant
        quadrants = zeros(2,2,2);
        % horizontal and vertical cuts
        cuts = [ceil(size(mask2D,1)/2), ceil(size(mask2D,2)/2)];
        helper = cell(2,2);
        helper{1,1} = mask2D(1:cuts(1),1:cuts(2));
        helper{1,2} = mask2D(1:cuts(1),cuts(2)+1:end);
        helper{2,1} = mask2D(cuts(1)+1:end,1:cuts(2));
        helper{2,2} = mask2D(cuts(1)+1:end,cuts(2)+1:end);

        quadrants(:,:,1) = cellfun(@(x) sum(sum(x)), helper);
        quadrants(:,:,2) = cellfun(@(x) sum(sum(x == 0)), helper);

        % project onto x- and y-axis
        projData = {sum(mask2D,1), sum(mask2D,2)};
        projData = cellfun(@(x) nnz(find(x == 0)), projData);

        % check if pattern is horizontal or vertical
    %     ver = sum(abs(quadrants(:,1,:) - quadrants(:,2,:)),1);
    %     hor = sum(abs(quadrants(1,:,:) - quadrants(2,:,:)),2);

        % determine if upper or left part (corresponds to isUpper) was sampled
        if(projData(2) > projData(1))     %hor(1,1,1) > ver(1,1,1) && hor(1,1,2) > ver(1,1,2))
            % horizontal pattern
            if(sum(quadrants(1,:,1),2) > sum(quadrants(2,:,1),2) && sum(quadrants(1,:,2),2) <= sum(quadrants(2,:,2),2))
                isUpper = true;
                isLower = false;
                checkBorder = 'last';
            else
                isUpper = false;
                isLower = true;
                checkBorder = 'first';
            end
            pfDimEst = pfDimIn(1);
        else
            % vertical pattern
            if(sum(quadrants(:,1,1),1) > sum(quadrants(:,2,1),1) && sum(quadrants(:,1,2),1) <= sum(quadrants(:,2,2),1))
                isUpper = true;
                isLower = false;
                checkBorder = 'last';
            else
                isUpper = false;
                isLower = true;
                checkBorder = 'first';
            end
            pfDimEst = pfDimIn(2);
        end

        % determine partial fourier dimension
        % just to be absolutely sure
        bord2D = cell(1,2);
        bord2D{1} = -ones(1,size(mask2D,2));
        bord2D{2} = -ones(1,size(mask2D,1));
        for i=1:size(mask2D,1)
            tmp = find(mask2D(i,:), 1, checkBorder);
            if(~isempty(tmp))
                bord2D{2}(i) = tmp;
            end
        end
        for i=1:size(mask2D,2)
            tmp = find(mask2D(:,i), 1, checkBorder);
            if(~isempty(tmp))
                bord2D{1}(i) = tmp;
            end
        end
        [n_el, xcenter] = cellfun(@(x) hist(x, min(x(x > 0)):max(x)), bord2D, 'UniformOutput', false);
        helpdim = cellfun(@(x,y) x(y == 0), n_el, xcenter, 'UniformOutput', false);
        if(any(~cellfun(@isempty, helpdim)))
            [helpdim{cellfun(@isempty, helpdim)}] = deal(0);
            [~, helpdim] = max(cell2mat(helpdim));
            pfDim = pfDimIn(~ismember([1 2],helpdim));
            if(pfDim ~= pfDimEst) % TODO: check if vertical and horizontal pattern detection still fails for large acceleration
                imagine(mask2D);
                warning('partialFourier()::detectPFdimCS(): Check determination of partial fourier dimension');
            end
        else
            pfDim = pfDimEst;
        end

        % create outputs
        smplPtrnSparse = logical(mask);
        smplPtrn = logical(mask);
        subs = {':',':',':'};
        if(pfn == 0)
            % partial fourier value not initialized in meas file
            bpoint = cellfun(@(x,y) x(y == max(y)), xcenter, n_el, 'UniformOutput', false);
            bpoint = bpoint{pfDimIn == pfDim};

            if(isLower)
                tmp = [1 2];
                subs{pfDimIn(pfDimIn == pfDim)} = bpoint:size(mask2D,tmp(pfDimIn == pfDim));
                smplPtrn(subs{:}) = true;
            else
                subs{pfDimIn(pfDimIn == pfDim)} = 1:bpoint;
                smplPtrn(subs{:}) = true;
            end
        else
            if(isLower)
                tmp = [1 2];
                subs{pfDimIn(pfDimIn == pfDim)} = round((1-pfn)*size(mask2D,tmp(pfDimIn == pfDim))):size(mask2D,tmp(pfDimIn == pfDim));
                smplPtrn(subs{:}) = true;
            else
                tmp = [1 2];
                subs{pfDimIn(pfDimIn == pfDim)} = 1:round(pfn*size(mask2D,tmp(pfDimIn == pfDim)));
                smplPtrn(subs{:}) = true;
            end
        end   
    end
else
    % full sampling pattern

    % init some helper variables
    pfDim    = 0;               % partial Fourier reduction dimension
    isUpper  = false;
    isLower  = false;
    isPartialFourier = false(Ndim, 1);

    for d = 1:Ndim
        subs = cell(1,Ndim);
        [subs{1:Ndim}] = deal(1);
        subs{d} = 1:sz(d);

        tmp = mask(subs{:});

        if(length(unique(tmp)) > 1)
            % partial fourier dimension found
            pfDim = d;
            isPartialFourier(d) = true;

            if(mask(1,1,1) == 1)
                isUpper = true;
            else
                isLower = true;
            end
            break;
        end
    end

    switch nnz(isPartialFourier)    % |N|umber of |N|on |Z|ero elements
        case 0
            error( 'pocs:NoPfDim', 'No partial Fourier dimension found.' )
        case 1
%             fprintf( 'Found partial Fourier along array dimension %d\n', pfDim + ~wasAddedCoilDim )
        otherwise
            error( 'pocs:TooManyPfDims', 'Partial Fourier only allowed in 1 dimension!' )
    end
    
    smplPtrn = logical(mask);
    smplPtrnSparse = logical(mask);
end
end


function [applDim, idx, idxBorder] = getCalibrationSize(ksp, mask, pfDim)

sz = size(mask);
Ndim = ndims( mask );
% applDim: mask applied along dimensions (0: not applied, 1: 1st dimension (pfDim),
% 2: 2nd dimension (loopDim)), structure: y-x-z
applDim = zeros(1,Ndim);
applDim(pfDim) = 1;

% determine the central line
tmp = squeeze( sum(abs(ksp),1) );
for d = 1:Ndim
   if d ~= pfDim
       tmp = max(tmp,[],d);     % keep only the maximum of non-partial data points
   end
end
[ ~, centreLine] = max( tmp(:) );


if(size(mask,3) > 1)
    % 3D
    
    % rotate mask to reduce testing complexity
    if(pfDim == 3)
        % z direction
        permRule = [ 3 1 2 ];
%         shiftRule = { sz(3):-1:1, ':' , ':' };
    elseif(pfDim == 2)
        % x direction
        permRule = [ 2 1 3 ]; 
%         shiftRule = { sz(2):-1:1, ':' , ':' };
    else
        % y direction
        permRule = [ 1 2 3 ];
%         shiftRule = { ':' , ':' , ':' };
    end
    
    helpermask = permute(mask,permRule); % it is ensured that upper part of mask is fully sampled
%     helpermask = helpermask(shiftRule{:}); % to ensure acquired kspace is lower part
    testmask = helpermask(:,:,1);
    
    idxBorder = -ones(1,size(testmask,2));
    for col=1:size(testmask,2)
        % get row index of border to fully sampled region
        tmp = find(testmask(:,col),1,'last');
        if(tmp)
            idxBorder(1,col) = tmp;
        end
    end
    
    if(any(idxBorder == -1))
        error('partialFourier(): unsampled column not possible for partial Fourier');
    end
    
    if(length(unique(idxBorder)) == 1)
        % flat border
        if(pfDim == 1 || pfDim == 2)
            % case 1 and 6
            applDim(3) = 2;
        else
            % case 3
            applDim(2) = 2;
        end
        
        % recalculate idxBorder in case of flat border 
        idx_helper = {':',':',':'};
        idx_helper{applDim == 0} = 1;
        testmask = squeeze(helpermask(idx_helper{:}));
        idxBorder = -ones(1,sz(applDim == 2));
        for col=1:sz(applDim == 2)
            % get row index of border to fully sampled region
            tmp = find(testmask(:,col),1,'last');
            if(tmp)
                idxBorder(1,col) = tmp;
            end
        end
    else
        if(pfDim == 3 || pfDim == 2)
            % case 4 and 5
            applDim(1) = 2;
        else
            % case 2
            applDim(2) = 2;
        end
    end
    
    
else
    % 2D
    loopDim = [1 2];
    loopDim = loopDim(~ismember(loopDim,pfDim));
    applDim(loopDim) = 2;
    
    if(loopDim == 1)
        helpermask = mask';
    else
        helpermask = mask;
    end
    
    idxBorder = -ones(1,size(helpermask,2));
    for col=1:size(helpermask,2)
        % get row index of border to fully sampled region
        tmp = find(helpermask(:,col),1,'last');
        if(tmp)
            idxBorder(1,col) = tmp;
        end
    end
    
    if(any(idxBorder == -1))
        error('partialFourier(): unsampled column not possible for partial Fourier');
    end
end  

halfLine = ceil(sz(applDim == 1)/2);
if(~ismember(halfLine,centreLine - 10:centreLine + 10))
    centreLine = halfLine;
end
   
idx = cell(1,sz(applDim == 2));
for col=1:sz(applDim == 2)
    dist = abs(idxBorder(1,col) - centreLine);
    idx{1,col} = centreLine - dist:centreLine + dist;
end
    

    


% sx = 2;
% sy = 2;
% 
% xflag = false;
% yflag = false;
% if(size(mask,3) > 1)
%     zflag = false;
%     sz = 2;
% else
%     zflag = true;
%     sz = 1;
% end
% 
% while(~(xflag && yflag && zflag))
% 
%     if(~xflag)
%         mask_helper = crop(mask,[sx+1,sy,sz]);
%         if(all(mask_helper))
%             sx = sx + 1;
%         else
%             xflag = true;
%         end
%     end
% 
%     if(~yflag)
%         mask_helper = crop(mask,[sx,sy+1,sz]);
%         if(all(mask_helper))
%             sy = sy + 1;
%         else
%             yflag = true;
%         end
%     end
%     
%     if(~zflag)
%         mask_helper = crop(mask,[sx,sy,sz+1]);
%         if(all(mask_helper))
%             sz = sz + 1;
%         else
%             zflag = true;
%         end
%     end
% 
%     if(sx == size(mask,1))
%         xflag = true;
%     end
%     if(sy == size(mask,2))
%         yflag = true;
%     end
%     if(sz == size(mask,3))
%         zflag = true;
%     end
% end
% 
% if(size(mask,3) > 1)
%     calibSize = [sx,sy,sz];
% else
%     calibSize = [sx,sy];
% end
% 
% m = size(mask);
% centers = floor(m/2);
% centers(pfDim) = centreLine;
% idx = cell(1,length(calibSize));
% for n=1:length(calibSize)
%     idx{n} = centers(n)+ceil(-calibSize(n)/2) : centers(n)+ceil(calibSize(n)/2);
% end

end
