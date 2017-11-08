%% demo code
% reconstruction of retrospective subsampled Shepp-Logan phantom

% set paths
addpath(genpath(pwd));
addpath(genpath(['..',filesep,'CS_LAB']));

%% ------------------------------------------
% 1.) load data: phantom (2 slices, 1 channel)
% -------------------------------------------
dImg = phantom('Modified Shepp-Logan',256);
dKSpace = fftnshift(dImg,1:2);


%% --------------------------
% 2.) create subsampling mask
% ---------------------------
if(ispc), sExt = '.exe'; else sExt = ''; end
system(['.',filesep,'sampling',filesep,'Subsample',sExt,sprintf(' %d %d %d %d', 256, 1, 4, 1)]);
iMask = readSamplingMaskFromFile(['.',filesep,'samplingPattern.txt']);


%% -------------------------
% 3.) apply subsampling mask
% --------------------------
% k-Space dimensions
% cell array: NSlice x NChannels x NRepetitions x NAverages
%             each cell: 2D: yPhase x xFreq x (nTime/nPhases)
%                        3D: yPhase x xFreq x zPhase
%                        4D: yPhase x xFreq x zPhase x nTime/nPhases
dKSpaceSub = {dKSpace .* repmat(iMask,[1 256])};


%% --------------------
% 4.) reconstruct image
% ---------------------
% set some recon parameters
para.cstype = 'FOCUSS';
para.transformation = 'fft';
para.lambda = 1e-12;
para.lambdaTV = 1e-6;
para.measPara.dim = [size(dKSpaceSub{1}), 1, 1, 1];
para.measPara.LCall = ones(1,4);
para.measPara.dimension = '2D';
para.espresso.state = false;
para.espresso.direction = 'off';
para.espresso.pfn = 1;
para.flagOversampling = logical([0 0 0]);
para.postproc.turnImage = false;
para.postproc.FreqOversamplingCorr = false;
para.prop.flagPlot = false;

% CS reconstruction
dCSRecon = CS_reconstruction(dKSpaceSub, para);

% zero-padded reconstruction
dZeropadded = abs(ifftnshift(dKSpaceSub{1}));


%% ------------------
% 5.) compare results
% -------------------
figure;
subplot(1,3,1)
imagesc(dImg); colormap('gray');
axis 'equal'; set(gca,'XTick',[],'YTick',[]); axis 'off';
title('original')
subplot(1,3,2)
imagesc(dZeropadded); colormap('gray');
axis 'equal'; set(gca,'XTick',[],'YTick',[]); axis 'off';
title('zero-padded recon')
subplot(1,3,3)
imagesc(dCSRecon); colormap('gray');
axis 'equal'; set(gca,'XTick',[],'YTick',[]); axis 'off';
title('CS recon')