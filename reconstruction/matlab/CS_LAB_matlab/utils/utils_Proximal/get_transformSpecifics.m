function [ TransformSpecifics ] = get_transformSpecifics( TransformName, reconDIM, Stages, FilterName, n1, n2, n3 )
%%
% Initializes filters and parameters for employed transform.

% M. Fischer, April 2016

%%

if ~exist('FilterName')
    TransformSpecifics.FilterName = 'db2';
else
    TransformSpecifics.FilterName = FilterName;
end;
TransformSpecifics.TransformName = TransformName;
TransformSpecifics.Stages = Stages;

if strcmp(reconDIM,'2D')
    image = zeros(n1,n2);
    [~, waveS] = wavedec2(image,Stages,FilterName);
    TransformSpecifics.waveS = waveS;
    
    switch TransformName
        case 'dwt'
            % see above
        case 'cplxdt_matlab'
            TransformSpecifics.CellAmount = TransformSpecifics.Stages+1;
            image = zeros(n1,n2);
            [dfFilter,rfFilter] = dtfilters('dtf2');
            TransformSpecifics.dfFilter = dfFilter;
            TransformSpecifics.rfFilter = rfFilter;
            wt_struct = dddtree2('cplxdt',image,3,dfFilter{1},dfFilter{2});
            TransformSpecifics.wt_struct = wt_struct;
        case 'cplxdt'
            intval = 2^TransformSpecifics.Stages;
            TransformSpecifics.extend_x = (mod(n1,intval)~=0)*intval-mod(n1,intval);
            TransformSpecifics.extend_y = (mod(n2,intval)~=0)*intval-mod(n2,intval);
            if TransformSpecifics.extend_x ~= 0 ||  TransformSpecifics.extend_y ~= 0
                TransformSpecifics.flagExtend = true;
            else
                TransformSpecifics.flagExtend = false;
            end;
            TransformSpecifics.CellAmount = TransformSpecifics.Stages+1;
            [Faf, Fsf] = FSfarras;
            [af, sf] = dualfilt1;
            TransformSpecifics.Faf = Faf;
            TransformSpecifics.af = af;
            TransformSpecifics.Fsf = Fsf;
            TransformSpecifics.sf = sf;
        case 'shearlet'
            if n1 <= n2
                TransformSpecifics.extend_x = n2-n1;
                TransformSpecifics.extend_y = 0;
            elseif n1 > n2
                TransformSpecifics.extend_y = n1-n2;
                TransformSpecifics.extend_x = 0;
            end;
            if TransformSpecifics.extend_x ~= 0 ||  TransformSpecifics.extend_y ~= 0
                TransformSpecifics.flagExtend = true;
            else
                TransformSpecifics.flagExtend = false;
            end;
            TransformSpecifics.lpfilt='maxflat'; % setup parameters for shearlet transform
            TransformSpecifics.CellAmount = 4+1; %TransformSpecifics.Stages+1;
            TransformSpecifics.shear_parameters.dcomp =[ 3  3  4  4]; % .dcomp(i) indicates there will be 2^dcomp(i) directions
            TransformSpecifics.shear_parameters.dsize =[32 32 16 16]; % .dsize(i) indicate the local directional filter will be dsize(i) by dsize(i)
            [~,TransformSpecifics.shear_f] = nsst_dec2(zeros(n1+TransformSpecifics.extend_x,n2+TransformSpecifics.extend_y), TransformSpecifics.shear_parameters, TransformSpecifics.lpfilt); % available versions: dec1, dec1e, dec2
        case 'shearlet_nontight'
            TransformSpecifics.shearletSystem = SLgetShearletSystem2D(0,n1,n2,Stages,ceil((1:Stages)/2),0); % add ,modulate2(dfilters('cd','d'),'c') to decrease amount of direction filters
        otherwise
            warning('No correct Transform chosen');
    end;
elseif strcmp(reconDIM,'3D')
    image = zeros(n1,n2,n3);
    TransformSpecifics.WDEC = wavedec3(image,Stages,FilterName);
    TransformSpecifics.WDEC.dec = [];
    
    switch TransformName
        case 'dwt'
            TransformSpecifics.CellAmount = 7*TransformSpecifics.WDEC.level+1;
            % see above
        case 'cplxdt'
            intval = 2^TransformSpecifics.Stages;
            TransformSpecifics.extend_x = (mod(n1,intval)~=0)*intval-mod(n1,intval);
            TransformSpecifics.extend_y = (mod(n2,intval)~=0)*intval-mod(n2,intval);
            TransformSpecifics.extend_z = (mod(n3,intval)~=0)*intval-mod(n3,intval);  
            if TransformSpecifics.extend_x ~= 0 ||  TransformSpecifics.extend_y ~= 0 ||  TransformSpecifics.extend_z ~= 0
                TransformSpecifics.flagExtend = true;
            else
                TransformSpecifics.flagExtend = false;
            end;
            TransformSpecifics.CellAmount = TransformSpecifics.Stages+1;
            [Faf, Fsf] = FSfarras;
            [af, sf] = dualfilt1;
            TransformSpecifics.Faf = Faf;
            TransformSpecifics.af = af;
            TransformSpecifics.Fsf = Fsf;
            TransformSpecifics.sf = sf;
        case 'shearlet' % see utils_ShearletTight\3Dshearlet_toolbox\DenoiseDemo\main.m
            % Due to nature of upsampling and downsampling for better denoising performance
            % dimension of data need to be divisible by 3*2^(number of decomposition level required).
            level=2; % choose level of decomposition % atm fixed since dBand has to be set also.
            % ! anything bigger than 2 takes ages.
            intval = 3*2^level;
            TransformSpecifics.extend_x = (mod(n1,intval)~=0)*intval-mod(n1,intval);
            TransformSpecifics.extend_y = (mod(n2,intval)~=0)*intval-mod(n2,intval);
            TransformSpecifics.extend_z = (mod(n3,intval)~=0)*intval-mod(n3,intval);  
            if TransformSpecifics.extend_x ~= 0 ||  TransformSpecifics.extend_y ~= 0 ||  TransformSpecifics.extend_z ~= 0
                TransformSpecifics.flagExtend = true;
            else
                TransformSpecifics.flagExtend = false;
            end;
            TransformSpecifics.CellAmount = 1+3*(4*4+2*2); % TransformSpecifics.Stages;
            % dataClass = 'single'; % 'single' or 'double'
            filterDilationType = '422'; %only two type '422' or '442'. currently only 422 supported
            filterType='meyer';%%only meyer type implemented
            % Cell to specify different directional band at as per level choosen.
            % In current implementation specifying second number of direction is ignored.
            dBand = {{[4 4]},{[4 4] [2 2]}};
            %dBand={{[8 8]}, ... % for level =1
            %    {[8 8],[6 6]}, ...  % for level =2
            %    {[8 8],[6 6],[4 4]}};   % for level =3
            %    {[8 8],[8 8],[4 4],[4 4]}}; % for level =4 
            filterSize=[12 12]; % [24 24 24 24]% choose filter size such that remainder is 0 when divided by number of direction in that level
            TransformSpecifics.shearingFilter = GetFilter(filterType,level,dBand,filterSize,filterDilationType ,'double'); %Build Windowing Filter for different Band
        case 'shearlet_nontight'
            if TransformSpecifics.Stages == 1
                TransformSpecifics.extend_z = max(0,13-n3); % < 13 will fail
            elseif TransformSpecifics.Stages == 2
                TransformSpecifics.extend_z = max(0,33-n3); % < 33 will fail
            else
                warning('Got enough memory? We will find out soon');
            end;
            if TransformSpecifics.extend_z ~= 0
                TransformSpecifics.flagExtend = true;
            else
                TransformSpecifics.flagExtend = false;
            end;
            % add modulate2(dfilters('cd','d'),'c') to decrease amount of direction filters
            % add floor((1:Stages)/2) to decrease amount of shearlets
            TransformSpecifics.shearletSystem = SLgetShearletSystem3D(0,n1,n2,n3+TransformSpecifics.extend_z,Stages,floor((1:Stages)/2),0,modulate2(dfilters('cd','d'),'c'));
        case 'surfacelet'
            % a faster alternative to 3D shearlets
            % see "Multidemensional Directional Filter Banks and Surfacelets or google "SurfBox".
        otherwise
            warning('No correct Transform chosen');
    end;
else
    error('ReconDim is not supported');
end;

end
   