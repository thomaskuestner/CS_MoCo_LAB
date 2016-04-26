%% Compressed Sensing
%%% optimization solver
% Zero-filling
% Zero: zero-filling and IFFT

% FOCUSS reconstruction
% FOCUSS: FOcal Underdetermined System Solver with conjugate gradient
%         solver

% sparseMRI
% M. Lustig, D. Donoho and JM Pauly "Sparse MRI: The Application
% of Compressed Sensing for rapid MR Imaging", Magnetic Resonance in
% Medicine, 2007
% sparseMRI: sparse nonlinear reconstruction with L1-norm minimization
%            with the help of a non-linear conjugate gradient

% SPIRiT and ESPIRiT reconstruction
% SPIRiT_*: M. Lustig and JM Pauly "SPIRiT: iTerative Self-consistent Parallel 
% Imaging Reconstruction from Arbitrary k-space", Magnetic Resonance in 
% Medicine, 2010
% SPIRiT_CG:   conjugate gradient SPIRiT reconstruction
% SPIRiT_POCS: projection onto convex sets SPIRiT reconstruction with
%              L1-sparsity
% ESPIRiT_*: Uecker et. al "ESPIRiT - An Eigenvalue Approach to Autocalibrating
% Parallel MRI: Where SENSE meets GRAPPA", Magnetic Resonance in Medicine,
% 2013 DOI 10.1002/mrm.24751
% ESPIRiT_CG:  conjugate gradient ESPIRiT reconstruction in k-space
% ESPIRiT_L1:  conjugate gradient ESPIRiT reconstruction with L1-wavelet
%              thresholding in image domain

% BART 
% BART: Berkeley Advanced Reconstruction Toolbox  

% L1-Magic
% L1_Magic_*: reconstruction via l1-magic package from Emmanuel Candes and 
% Justin Romberg, Caltech: "l1-magic: Recovery of Sparse Signals via 
% Convex Programming", October 2005
% L1_Magic_TV:        TV minimization
% L1_Magic_L1:        L1 minimization with quadratic constraints
% L1_Magic_TVDantzig: Dantzig TV minimization 
% L1_Magic_L1Dantzig: L1 minimization with bounded residual correlation

% Proximal averages 
% uses wavelet_mat, in brackets reconstruction dimensionality
% and number space
% synthesis model:
% *_proxA (Proximal Average): Average of proximal operators. Based on Y. Yu 
% "Better Approximation and Faster Algorithm using the Proximal Average", 
% 2013.
% FCSA_* : Implementation according to Huang et al. "Efficient MR image 
% reconstruction for compressed MR imaging", 2011.
% FCSA_WaTMRI: implementation according to Huang et al. "Exploiting the 
%              wavelet structure in compressed sensing MRI", MRI 2014      (2D real)
% FCSA_SLEP:   WaTMRI modifications for group sparsity using SLEP (Liu et
%              al., "SLEP: Sparse Learning with Efficient Projections",
%              2009) to solve TSGLasso                                     (2D real)
% FCSA_proxA:  FCSA with proximal averages                                 (2D real/complex, 3D real/complex)
% BFCSA_proxA: FCSA with an outer Bregman Iteration based on BOS method by 
%              Zhang et al. "Bregmanized nonlocal regularization for 
%              deconvolution and sparse reconstruction", 2012              (2D real/complex, 3D_sliced real/complex)
% ADMM_proxA:  ADMM with proximal average according to Afonso et al. "Fast 
%              Image Recovery Using Variable Splitting and Constrained 
%              Optimization", 2010                                         (2D real/complex, 3D_sliced real/complex) 
% SB_proxA:    Split-Bregman according to Goldstein et al. "The Split  
%              Bregman Method for L1-regularized Problems", 2009           (2D real/complex, 3D_sliced real/complex)
% analysis model:
% A_PFISTA:    projected FISTA according to Liu et al., "Projected 
%              Iterative Soft-Thresholding Algorithms for linear inverse 
%              problems", 2015                                             (2D real/complex, 3D real/complex)
% A_TDIHT:     transform domain IHT according to Giryes, "A greedy 
%              algorithm for the analysis transform domain", 2016          (2D real/complex, 3D real/complex)
% A_ADMM_L1:   analysis ADMM with l1 penalty                               (2D real/complex)
% A_ADMM_L0:   analysis ADMM with l0 penalty                               (2D real/complex)
% A_ADMM_SCAD: analysis ADMM with Smoothly Clipped Absolute Deviation
%              penalty according to Fan and Li, "Variable selection via 
%              nonconcave penalized likelihood and its oracle properties",
%              2001                                                        (2D real/complex, 3D real/complex)
% A_ADMM_MCP:  analysis ADMM with Minimax Concave Penalty according to
%              Zhang, "Nearly unbiased variable selection under minimax 
%              concave penalty", 2010                                      (2D real/complex, 3D real/complex)
% A_ADMM_ATAN: analysis ADMM with arctangent penalty according to Selesnick
%              and Bayram, "Sparse signal estimation by maximally sparse 
%              convex optimization", 2014                                  (2D real/complex, 3D real/complex)
% A_ADMM_PSHRINK: analysis ADMM with p-shrinkage penalty according to 
%                 Woodworth and Chartrand, "Compressed sensing recovery via
%                 nonconvex shrinkage penalties", 2015                     (2D real/complex, 3D real/complex)

cstype = 'FOCUSS';


%%% sparsifying transformation
% fft:          Fourier transformation
% pca:          Principal component analysis
% dct:          Discrete cosine transformation
% wavelet_mat:  Wavelet (via Mathworks Wavelet toolbox)
% wavelet_lab:  Wavelet (via WaveLab850) - ATTENTION: just 2D possible!
% mellin:       Fourier-Mellin transformation
% gabor:        Gabor transformation (TODO)
% surfacelet:   Surfacelet transformation
% bandlet:      Bandelet transformation (TODO)
% curvelet:     Contourlet transformation
% dict:         Dictionary learning (TODO)
transformation = 'fft';


%% set algorithm parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% k-t FOCUSS parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outer loops, default:
% FOCUSS: 2
% sparseMRI: 1
% proximal averages: 500
iNOUTER = 2;
% inner loops (CG), default:
% ktFOCUSS: 20
% SPIRiT: 30
% ESPIRiT: 5
% sparseMRI: 45
% RecPF: 50
% proximal averages: 20
iNINNER = 20;
% iteration limit for NLTV (default: 120)
itrNLTV = 120;

% use l_p norm (default: 0.5)
p = 0.5;
% sparsifying transformation (default:5, max: 75)
lambda = 1e-4; %0.01; 5
% emphasize calibration consistency (default: 50, max: 350)
lambdaCalib = 0; %0.05; 50 
% emphasize total variation (TV) (default: 0.001)
lambdaTV = 0;
% emphasize ESPReSSo conjugate similarity
lambdaESPReSSo = 0; %0.02; 10
% emphasize motion correction
lambdaMC = 0;
% group sparsity (default: 0.005)
lambdaGroup = 0.005;
% NLTV regularizer (lambda and sigma weighting)
lambdaNLTV = 0.5; % default: 0.5
lambdaNLTV_h = 2.6; % may fail if chosen too small, default: 0.001 (WaTMRI), 2.6 (ADMM/SB)
% step size for ADMM, SB/tradeoff for WaTMRI
mue = 0.2;
% relative regularizer weights for averaging in CSA 
regularizerWeights = [1 1 1 1]; % [Wavelet, TV, Group/Tree-Wavelet, NLTV/NLM]

% noise level (default: 1e-6, 1e-9 (for L1_TV))
epsilon = 1e-6;


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Proximal Average %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% regularizer flags
% total variation
flags.flagTV = true;
% isotropic total variation
flags.flagTV_iso = true;
% wavelet 
flags.flagWave = true;
% group sparsity
flags.flagGroup = true;
% NLTV
flags.flagNLTV = false;
% SB NLTV
flags.flagSBNLTV = true;

% reconstruction flags
% real, complex processing
flags.flagRealAndImag = true; % separate into real and imaginary
flags.flagRealOnly = false;

% set reconstruction dimensionality for proximal averages
% 2D | 3D (just available for 3D images) - will be set automatically
reconDIM = '2D';

% apply sense mask for multichannel reconstruction
lSense = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FOCUSS parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scaling parameter for weighting matrix:
% 'none' | 'self' | 'energy'
sScaling = 'energy';

% opt problem (OLD -> throw out!!!)
opt_problem = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GRAPPA / SPIRiT parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calibration kernel size
% 2D: k_y - k_x
% 3D: k_y - k_z - k_x
kernelSize = [7,7,9];

% indicate whether kernelImg should be zero padded to
% size(img)+kernelSize-1 (true) or not (false)
% +: suppresses aliasing, -: more computational effort
flagZeropadding = false;

% Tykhonov regularization in the calibration (default: 0.01)
calibTyk = 0.01;  

% Tykhonov regularization in the reconstruction (SPIRiT/ESPIRiT) (default: 1e-5)
reconTyk = 1e-5;  

% set a fixed calibration size OR leave blank for automatic calculation
% 2D: k_y - k_x (k_y: middle = floor(nPha/2); y = middle-4+1:middle+4; 
%                k_x: complete length)
% 3D: k_y - k_z - k_x
% calibSize = [11, 11, nFreq];
calibSize = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ESPIRiT parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% either use ESPIRiT toolbox or MatLab implementation
flagToolbox = true;

% number of splitting iterations for CS part
iNIterSplit = 15;    

% splitting parameter (default: 0.4 or 0.5)
splitWeight = 0.5;

% threshold of eigenvectors in k-space (default: 0.02)
eigThresh_k = 0.002; 

% threshold of eigenvectors in image space (default: 0.9)
eigThresh_im = 1e-4; 

% weight the sensitivity maps with n maps of eigenvalues (default: 2)
n_maps = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NLCG/sparseMRI %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% l1 smoothing parameter
l1Smooth = 1e-15;

% backtracking line search iteration limit
lineSearchItnlim = 150;

% backtracking line search step sizes alpha and beta
lineSearchAlpha = 0.01;
lineSearchBeta = 0.6;

% backtracking line search start point
lineSearchT0 = 1;

% NLCG stopping condition
gradToll = 1e-9;

%%%%%%%%%%%%%%%%%%
%%%% L1-Magic %%%%
%%%%%%%%%%%%%%%%%%
% log barrier algorithm terminates when the duality gap <= lbtol (default: 1e-6)
lbtol = 1e-6;

% increase of barrier constant at each iteration (default: 10)
tvmu = 10;

% tolerance for Conjugate Gradients (default: 1e-10)
cgtol = 1e-10;

% maximum number of iterations for Conjugate Gradients (default: 200)
cgmaxiter = 200;

% maximum number of iterations for primal dual algorithm (default: 50)
pdmaxiter = 50;


%%%%%%%%%%%%%%%%%%%%%%
%%%% oversampling %%%%
%%%%%%%%%%%%%%%%%%%%%%
% perform correction of x oversampling before reconstruction and y/z
% oversampling after reconstruction (due to backward compatibility if 
% drecksMDH is not available) (x-y-z/slice)
flagOversampling = logical([0 0 0]);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% cut-out window %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% cut out window
% no para:                @rectwin, @barthannwin, @bartlett, @bohmanwin, @parzenwin, @triang
% 'symmetric'|'periodic': @blackman, @blackmanharris, @hamming, @hann,
%                         @flattopwin, @nuttallwin, @lanczos, @welch, @cosine
% symmetry + para:        @chebwin (60), @tukeywin (0.25), @gausswin (0.3), @kaiser (4), @planckTaper (0.25)
FFTwindow.type = @rectwin;
FFTwindow.windowOpt = {'symmetric', 4};

%%%%%%%%%%%%%%%%%%
%%%% ESPReSSo %%%%
%%%%%%%%%%%%%%%%%%
% reconstruction type: inside of CG = 1, just outside of CG = 2
espresso.reconType = 1;
% POCS iterations of ESPReSSo recon
espresso.iter = 30;
% espresso norm: 1 = L1, 2 = L2
espresso.norm = 2;
% espresso constraint: 1 = POCS, 2 = Im(rho)
espresso.constraint = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Motion Correction %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path containing elastix and transformix (no white spaces in path allowed)
mc.sElastixPath = '';
% SURF landmark features
mc.lSURF = false;
% use SENSEmap in recon for correct RSS->channel mapping
mc.lSenseMC = true;
% parameter files
% normal: single B-spline registration
% concatenated: concatenated rigid, affine and B-Spline registration 
mc.sElastixParams = 'normal';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% postprocessing image reconstruction %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rss: root-sum-of-squares
% SENSE: SENSE based channel addition
postproc.type = 'rss';
% set tolerance value for rank estimation of SENSE
postproc.tol = 0.004;
% correct for anisotropic pixel resolution (y-x-z/slices-t)
postproc.aniso = [];
% anisotropic and (rotation) interpolation method ('nearest', 'bilinear')
postproc.interp = 'bilinear';
% switch image turning on/off
postproc.turnImage = true;
% correct frequency oversampling
postproc.FreqOversamplingCorr = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% sparsifying transformation parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trafo.trafoType = transformation;
if(nargin > 1 && exist('para', 'var') && isfield(para, 'transformation')), trafo.trafoType = para.transformation; end;

%%%% Transform dimensions %%%%
% transform dimensions during CS recon
% starting input:
% 2D: k_y-k_x or t-k_y-k_x
% 3D: k_y-k_z-k_x
% 4D: t-k_y-k_z-k_x
% 1st line: indicates wether the input image should be fourier transformed
% DURING CS reconstruction from k_A -> A (if value is "1") or not ("0").
% For applying an inverse fourier transformation (k_A -> A) BEFORE the CS
% reconstruction, set the value to "2"
% 2nd line: indicates (with "1") which domains should be transformed into the new
% basis (will be automatically flipped to the preceding dimensions)
% use logical arrays assuming the following order: [y x z t]
% sparse dimensions: y, z, t
% example a): 
% [1 0 1 1; 1 1 0 0] 
% -> fft on t, y, z during compBasis (1st line);
% sparsifying transformation (2nd line) on y (y_cart -> y_newBase) and 
% x (kx_cart -> kx_newBase), all other dimensions will be left untouched, here t ands z
% example b):
% [0 1 1 0; 1 0 1 1]
% -> fft on x and z during compBasis (1st line) (forward: kx -> x, backward x -> kx);
% sparsifying transformation (2nd line) on y (ky_cart -> ky_newBase), 
% z (z_cart -> z_newBase) and t (t_cart -> t_newBase)
trafo.trafodim = [1 2 1 0; 1 0 1 0];

% set dimensions to scramble for fourier trafo
% either the same as the fftdim or 
trafo.scrambledim = trafo.trafodim(1,:);


% specify whether the image should be rescaled (false) or zero-padded (true) 
% if necessary (wavelet_lab, pca, surfacelet and kernelMult for flagZeropadding=true)
trafo.zeroPad = true;
% interpolation method for image rescaling ('nearest', 'bilinear',
% 'bicubic', 'box', 'triangle', 'cubic', 'lanczos2', 'lanczos3')
trafo.rescaleInterp = 'bilinear';

switch trafo.trafoType
    case 'fft'
        trafo.windowing = false;
        % no para:                @rectwin, @barthannwin, @bartlett, @bohmanwin, @parzenwin, @triang
        % 'symmetric'|'periodic': @blackman, @blackmanharris, @hamming, @hann,
        %                         @flattopwin, @nuttallwin, @lanczos, @welch, @cosine
        % symmetry + para:        @chebwin (60), @tukeywin (0.25), @gausswin (0.3), @kaiser (4), @planckTaper (0.25)
        trafo.windowType = @rectwin;
        trafo.windowOpt = {'symmetric', 4};
                
    case 'wavelet_mat'
        % wavelet filter and corresponding size
        % 'haar':   Haar (='db1')
        % 'db':     Daubechies (1:45)
        % 'coif':   Coiflets (1:5)
        % 'sym':    Symlets (2:45)
        % 'dmey':   Discrete Meyer 
        % 'bior':   Biorthogonal (1.1, 1.3, 1.5, 2.2, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.7, 3.9, 4.4, 5.5, 6.8)
        % 'rbio':   Reverse Biorthogonal (1.1, 1.3, 1.5, 2.2, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.7, 3.9, 4.4, 5.5, 6.8)
        trafo.waveletFilter = 'db';
        trafo.waveletFilterSize = 2;
        trafo.waveletFilter_l12 = 'db'; % for group sparsity
        trafo.waveletFilterSize_l12 = 1;
        % extension mode 
        % 'zpd':            zero padding
        % 'sym'|'symh':     Symmetric-padding (half-point)
        % 'symw':           Symmetric-padding (whole-point)
        % 'asym'|'asymh':   Antisymmetric-padding (half-point)
        % 'asymw':          Antisymmetric-padding (whole-point)
        % 'spd'|'sp1':      Smooth-padding of order 1
        % 'sp0':            Smooth-padding of order 0
        % 'ppd':            Periodic-padding
        trafo.extMode = 'zpd';
        % #decomposition stages
        trafo.waveletStages = 2;
        % apply soft thresholding 
        trafo.flagThresholding = false;
        % threshold value for soft threshold 
        % if empty: calculate value automatically!
        trafo.wavWeight = 0.0015;
        
    case 'wavelet_lab'
        % qmf filter type:
        % 'Haar', 'Beylkin', 'Coiflet', 'Daubechies', 'Symmlet', 'Vaidyanathan', 'Battle'
        % use suffix _TI for translation invariance, e.g.: 'Daubechies_TI'
        % TIflag: 0 orthogonal, 1 translation invariant
        trafo.waveletFilter = 'Daubechies';
        % qmf filter size
        % wavelet has p vanishing moments
        % more vanishing moments:
        % +) complex functions can be represented with a sparser set of wavelet coefficients
        % -) higher computational burden
        % depending on the filter type the filter size sets the amount of vanishing
        % moments p (e.g. Daubechies: FilterSize=a => p=a/2 vanishing moments)
        % valid filter sizes:
        % Haar, Beylkin, Vaidyanathan: arbitrary (>0)
        % Coiflet: 1:5
        % Daubechies: 4:2:20
        % Symmlet: 4:10
        % Battle: 1,3,5
        trafo.waveletFilterSize = 4;
        % coarsest level-l of decomposition := #filter stages
        % size of coarsest image (approximation) in level-l: diaSize/2^l
        % => input to wavelab850: L = log2(diaSize/2^waveletStages)
        % wavelet image will have 2^N x 2^N pixels in the coarsest level
        % the smaller N, the more filter stages (#decompositions)
        trafo.waveletStages = 4;        
        % soft-thresholding value (default: SPIRiT: 0.0015, ESPIRiT (PC/Unix): 1e-9/0.0015)
        trafo.wavWeight = 0.00015;
        
        % apply soft thresholding (always true for SPIRiT/ESPIRiT)
        trafo.flagThresholding = true;
        
    case 'mellin'
        % distance scaling factor for AFMT
        trafo.sigma = 0.5; 
        % value used to extrapolate pixels between circumscribed circle and 
        % the border of the image in cartesian coordinates
        trafo.extrapVal = 0;
        % size of the AFMT in multiples of the original image dimensions
        % (1st dimension, 2nd dimension)
        trafo.trafoSize = [2, 2];
        % data interpolation method
        % 'nearest':    nearest neighbour interpolation
        % 'linear':     linear interpolation
        % 'cubic':      cubic interpolation -> not supported
        % 'spline':     spline interpolation -> not supported
        trafo.interp = 'linear';
        % backward transformation type
        % 'fast':       interpolate in r-phi plane (fast but inaccurate)
        % 'accurate':   interpolate in cartesian plane (accurate but slow)
        trafo.bTrafoType = 'fast';
        
    case 'curvelet'
        % number of radial scales
        trafo.nbscales 	    = 3;    
        % use wavelet for finest scale
        trafo.allCurvelets  = 0;  
        % number of angular scales (has to be a multiple of 4)
        trafo.nbdstz_coarse = 8;     

    case 'surfacelet'
        % type of multiscale pyramid, corresponding to different levels
        % of redundancy
        % 1:    ~ 6.4
        % 1.5:  ~ 4.0
        % 2:    ~ 3.4
        trafo.Pyr_mode = 2;
        % filter name for the hourglass filter bank
        % 'ritf': rotational-invariant tight frame
        trafo.HGfname = 'ritf';
        % order of checkerboard filters (default: 12)
        trafo.bo = 12;
        % size of the mapping kernel. This controls the quality of the
        % hourglass filters (default: 15)
        trafo.msize = 15;
        % beta value of the Kaiser window (default: 2.5)
        trafo.beta = 2.5;
        % parameter lambda used in the hourglass filter design (default: 4)
        trafo.lambda = 4;
        % pad all matrices to the same maximal size (if false: forward and
        % backword transformation need to be calculated for correct zero-padding and cropping)
        trafo.padOnSameSize = true;
        % #decomposition levels L (default: 2)
        trafo.decompLevels = 2;
        % further subband decomposition at each level l
        % diagonal element (i,i) must be -1: hour-glass filter along
        % i-direction
        % off-diagonal element (i,j): hour-glass filter along i-direction
        % is decomposed into 2^x subbands along j-direction
        trafo.Lev_array = cell(1,trafo.decompLevels);
        [trafo.Lev_array{:}] = deal(diag(-1*ones(1,3))); % 3 := 3D -> other dimensionality will be corrected in checkConsExist
        trafo.Lev_array{end} = ones(3);
        trafo.Lev_array{end}(logical(eye(size(trafo.Lev_array{end})))) = -1;
        % or if Lev_array contains scalar values (1 := dual-tree real
        % wavelet decomposition, 2 := dual-tree complex wavelet
        % decomposition)
               
    case 'bandelet'
        
    case 'gabor'
        
    case 'dict'
end


%% set calculation precision
if((exist('measPara','var') && isempty(measPara)) || ~exist('measPara','var') || (exist('measPara','var') && isstruct(measPara) && ~isfield(measPara,'precision')))
    measPara.precision = 'single';
end


%% set optimization for FFT
prop.flagOptim = true;
% planner method:
% 'measure', 'patient', 'exhaustive'
fft_planner_method = 'exhaustive';

%% check automatically for non-accelerated image
lAutoCheck = true;