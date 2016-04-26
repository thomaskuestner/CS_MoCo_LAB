function y = denoising_dtdwt(x)
% Local Adaptive Image Denoising Algorithm
% Usage :
%        y = denoising_dtdwt(x)
% INPUT :
%        x - a noisy image
% OUTPUT :
%        y - the corresponding denoised image

% Set the windowsize and the corresponding filter
windowsize  = 7;
windowfilt = ones(1,windowsize)/windowsize;

% Number of Stages
J = 6;
I=sqrt(-1);

% symmetric extension
L = length(x); % length of the original image.
N = L+2^J;     % length after extension.
x = symextend(x,2^(J-1));

load nor_dualtree    % run normaliz_coefcalc_dual_tree to generate this mat file.

% Forward dual-tree DWT
% Either FSfarras or AntonB function can be used to compute the stage 1 filters  
%[Faf, Fsf] = FSfarras;
[Faf, Fsf] = AntonB;
[af, sf] = dualfilt1;
W = cplxdual2D(x, J, Faf, af);
W = normcoef(W,J,nor);


% Noise variance estimation using robust median estimator..
tmp = W{1}{1}{1}{1};
Nsig = median(abs(tmp(:)))/0.6745;

for scale = 1:J-1
    for dir = 1:2
        for dir1 = 1:3
            
            % Noisy complex coefficients
            %Real part
            Y_coef_real = W{scale}{1}{dir}{dir1};
            % imaginary part
            Y_coef_imag = W{scale}{2}{dir}{dir1};
            % The corresponding noisy parent coefficients
            %Real part
            Y_parent_real = W{scale+1}{1}{dir}{dir1};
            % imaginary part
            Y_parent_imag = W{scale+1}{2}{dir}{dir1};
            % Extend noisy parent matrix to make the matrix size the same as the coefficient matrix.
            Y_parent_real  = expand(Y_parent_real);
            Y_parent_imag   = expand(Y_parent_imag);
            
            % Signal variance estimation
            Wsig = conv2(windowfilt,windowfilt,(Y_coef_real).^2,'same');
            Ssig = sqrt(max(Wsig-Nsig.^2,eps));
            
            % Threshold value estimation
            T = sqrt(3)*Nsig^2./Ssig;
            
            % Bivariate Shrinkage
            Y_coef = Y_coef_real+I*Y_coef_imag;
            Y_parent = Y_parent_real + I*Y_parent_imag;
            Y_coef = bishrink(Y_coef,Y_parent,T);
            W{scale}{1}{dir}{dir1} = real(Y_coef);
            W{scale}{2}{dir}{dir1} = imag(Y_coef);
            
        end
    end
end

% Inverse Transform
W = unnormcoef(W,J,nor);
y = icplxdual2D(W, J, Fsf, sf);

% Extract the image
ind = 2^(J-1)+1:2^(J-1)+L;
y = y(ind,ind);