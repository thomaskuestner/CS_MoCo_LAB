function fCreateLoopImages(dIn,SImg)

if(nargin < 2)
   lDF = false;
else
   lDF = true;
end

sTmpDir = 'C:\tmp';

lSave = true;
% show image
iSlice = 34;
dScale = [0 5.0E-8];
mcPath = '';
mcParaPath = ''; % 'elastix\para\normal';
sElastixPath = '';

%% always image
sFiles = dir([sTmpDir,filesep,'Bilder']);
iFiles = length(sFiles) - 2 + 1;

hfig = figure;
fullscreen = get(0,'ScreenSize');
set(hfig,'position',[fullscreen(1) fullscreen(2)+70 fullscreen(3) fullscreen(4)-150]);
set(hfig, 'Color', 'w');
iImg = uint8(dIn(:,:,iSlice,1)./(dScale(2) - dScale(1)).*255) + 1;
imagesc(iImg);
axis 'image';
colormap(gray(256));
axis 'off';
export_fig(hfig, [sTmpDir,filesep,'export',filesep,sprintf('img_%02d.tiff',iFiles)]);

if(lSave)
    save([sTmpDir,filesep,'Bilder',filesep,sprintf('img_%02d.mat',iFiles)],'dIn');
end
close(hfig);


%% deformation field
if(lDF)
    sFiles = dir([sTmpDir,filesep,'DF']);
    iFiles = length(sFiles) - 2 + 1;

    SDeform = fMCReconMain(dIn, SImg, sElastixPath, mcPath, mcParaPath, '', 0, true); % compose transformations

    hfig = fEvalRegistrationGray(iImg, iImg, SDeform(end).dFx(:,:,iSlice), SDeform(end).dFy(:,:,iSlice));
    export_fig(hfig, [sTmpDir,filesep,'export',filesep,sprintf('df_%02d.tiff',iFiles)]);

    if(lSave)
        save([sTmpDir,filesep,'DF',filesep,sprintf('df_%02d.mat',iFiles)],'SDeform');
    end
    close(hfig);
end


end
