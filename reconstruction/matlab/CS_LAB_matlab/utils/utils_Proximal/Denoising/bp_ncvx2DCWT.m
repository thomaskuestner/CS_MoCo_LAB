function [x] = bp_ncvx2DCWT(y,A,AH,J,lam,a,mu,Nit,pen)
% function [x] = bp_ncvx2DCWT(y,A,AH,J,lam,a,mu,Nit,pen)
%
% Convex denoising of 2D Image using non-convex tight frame regularization
% using the 2D dual tree complex wavelet transform
% Input:  
%     y - Input image
%     A - Forward transform (Undecimated Wavelet transform)
%     AH - Inverse transform
%     J - Number of scales
%     lam - Regularization parameter (vector)
%     a - Degree of non-convexity (a_i < 1/lam_i)
%     mu - Augmented Lagrangian parameter (mu > 1/r)
%     Nit - Number of iterations
%     pen - Regularizer to be used 
%             a. Logarithmic ('log')
%             b. Rational    ('rat')
%             c. Arctangent  ('atan')
%             d. L1          ('l1')
%     
% Output:
%     x - Denoised image
%         
% Ankit Parekh (ankit.parekh@nyu.edu), NYU School of Engineering
% Reference:
% Convex denoising using non-convex tight frame regularization
% Ankit Parekh and Ivan W. Selesnick
% IEEE Signal Process. Lett., 2015

x = zeros(size(y));
d = A(y);
k = 1/(1 + mu);
I = sqrt(-1);

for i = 1:Nit
    w = A(x);
    for j = 1:J
        for s1 = 1:2
            for s2 = 1:3
                C = (w{j}{1}{s1}{s2} + d{j}{1}{s1}{s2}) + I * (w{j}{2}{s1}{s2} + d{j}{2}{s1}{s2});
                C = thresh(C,lam/mu,a,pen);
                w{j}{1}{s1}{s2} = real(C);
                w{j}{2}{s1}{s2} = imag(C);
            end
        end
    end
                        
    for j = 1:J
        for s1 = 1:2
            for s2 = 1:3
                U{j}{1}{s1}{s2} = (w{j}{1}{s1}{s2} - d{j}{1}{s1}{s2});
                U{j}{2}{s1}{s2} = (w{j}{2}{s1}{s2} - d{j}{2}{s1}{s2}); 
            end
        end
    end
    U{J+1} = w{J+1};
    
    x = k * (y + mu * AH(U));
    u = A(x);
    for j = 1:J
        for s1 = 1:2
            for s2 = 1:3
                d{j}{1}{s1}{s2} = d{j}{1}{s1}{s2} - (w{j}{1}{s1}{s2} - u{j}{1}{s1}{s2});
                d{j}{2}{s1}{s2} = d{j}{2}{s1}{s2} - (w{j}{2}{s1}{s2} - u{j}{2}{s1}{s2});
            end
        end
    end
end

