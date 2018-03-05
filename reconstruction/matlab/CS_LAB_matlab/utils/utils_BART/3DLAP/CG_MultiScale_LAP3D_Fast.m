function u_est = CG_MultiScale_LAP3D_Fast(I1, I2, FilterSizes, PreFilt, MedFilt, uin)
% The function implements a fast multi-scale framework for the 3D version 
% of the LAP optical flow algorithm. Instead of downsampling the images, the 
% framework changes the size of the all-pass filters used in the LAP 
% algorithm. The filter basis used in the LAP algorithm spans the 
% derivatives of a Gaussian filter. 
%
% Function uses faster interpolation (but less accurate).
%
% Note that this implementation is for greyscale images only.
% 
% Input parameters:
%       I1 and I2   -> Input images of size M by N (Greyscale)
%       FilterSizes -> Controls the size of the filters (vector with
%                      values for R
%                      i.e. filter size is [2R + 1] by [2R + 1]
%       PreFilt     -> Optional parameter to decide whether to highpass
%                       filter the images. (Default = 1 (Yes), 0 = (No))
%       MedFilt     -> Optional parameter to decide whether to median
%                       filter the flow at fine scales. (Default = 1 (Yes), 0 = (No))
%
% Outputs:
%       u_est       -> Estimate of the optical flow 
%        
%
% Date: 08/07/2015, Author: Chris Gilliam
%
% modifications for MRI application:
% - different input parameter: Filter size instead of upper parameter R
% - initialization in single precision (except filter kernel)
% 26.08.2015, Verena Neumann/Thomas Kuestner

if nargin <= 3,
    PreFilt = 1;
    MedFilt = 1;
elseif nargin <= 4,
    MedFilt = 1;
end

% Obtain the dimensions of the images:
[M, N, P] = size(I1);

% Initialise filter functions:
funs = Filter_Functions;

% Local Gaussian filter (used in high pass filtering)
h = exp(-(-2:2).^2/2);
h = h./sum(h(:));

% Initial optical flow estimate
if(nargin < 6)
    u_holder = repmat({zeros(M,N,P,'single')}, 1, 3);
else
    u_holder = uin;
end

% define half support of filters (i.e. R)
amp_array = FilterSizes;

% Initialise local counter:
num_level = 0;

% Local I1 variable:
if PreFilt == 1,
    im1 = I1 - funs.Filter_General(I1, h); % High pass filtering using gaussian filter h
%       im1 = funs.Filter_Laplacian(im1);  % high-pass using generalised Laplacian filter
else
    im1 = I1;
end

% Start estimating the optical flow
for l = 1:length(amp_array),
    num_level = num_level + 1;
    disp(['Level ', int2str(num_level), '/', int2str(size(amp_array,2))]);
    
    % Define filter parameter R at each iteration
    amp_size = amp_array(l);
   
    % Load Filter Basis:
    Basis_Set = loadbasis(3,amp_size);
    
    tic;
    if l == 1,
        % No warping for first iteration:
        if PreFilt == 1,
            I2_shift = I2 - funs.Filter_General(I2, h); % High pass filtering using gaussian filter h
%             I2_shift = funs.Filter_Laplacian(I2);       % high-pass using generalised Laplacian filter
        else
            I2_shift = I2;
        end
    else
        % Warp I2 closer to I1 using current optical flow estimate
%         I2_shift = imshift_3D(I2,{-u_holder{1}, -u_holder{2}, -u_holder{3}}, 'shiftedlinear');
        I2_shift = ShiftedLinear_Interp_3D(-u_holder{1},-u_holder{2},-u_holder{3},I2);
        if PreFilt == 1,
            I2_shift = I2_shift - funs.Filter_General(I2_shift, h); % High pass filtering using gaussian filter h
%           I2_shift = funs.Filter_Laplacian(I2_shift); % high-pass using generalised Laplacian filter  
        end
    end
    
    disp(['Filter size = ', num2str(2*amp_size+1), ', Image Warping time = ', num2str(toc,3)]);
   
    % Using basis functions estimate optical flow on a local scale:
    tic;
    [uest_Orig, ~] = optiflowFilter3D(im1, I2_shift, Basis_Set);
    disp(['Filter size = ', num2str(2*amp_size+1), ', LAP Algorithm time = ', num2str(toc,3)]);

    % Clean optical flow
    % Stage 1: Remove nan's in the optical flow using inpainting:
    if sum(isnan(uest_Orig{1}(:))) >= M*N*P,
        error('All NaN. Suggest reduce Level_Num');
    end
    tic;
    [uest_Clean, ~] = cleanOF3D(uest_Orig);

    % Stage 2: Remove flow elements that corresponding to large warping
    % errors. 
    R = round(2*amp_size);
    k1 = -R:R;
    uest_Clean = cellfun(@funs.Filter_Gauss, uest_Clean, {R,R,R}, {k1,k1,k1}, 'uni',false);
        
    disp(['Filter size = ', num2str(2*amp_size+1), ', Post-Processing time = ', num2str(toc,3)]);
    
    % Rescale optical flow and add to estimate
    u_holder = cellfun(@plus,u_holder,uest_Clean,'uniformoutput',false);

    % Refinement of OF at highest level based on errors in shifted 
    %           images 
    if amp_size <= 2 && MedFilt == 1,
        u_holder = Cleaning_Procedure(u_holder);
    end

end 
        
u_est = u_holder;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Embedded functions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_out = Cleaning_Procedure(u_in)
% function cleans the estimate of the optical flow. The cleaning process
% comprises a robust smoothing stage using two Median filters (of different
% sizes)
%

% Define size of median filters
B1 = 11;
B2 = 3;
        
% Two part median filtering:
% fine scale
u_out = cellfun(@medfilt3, u_in, {B2,B2,B2}, {'symmetric','symmetric','symmetric'},'uni',false);
% coarse scale
u_out = cellfun(@medfilt3, u_out, {B1,B1,B1}, {'symmetric','symmetric','symmetric'},'uni',false);
end

