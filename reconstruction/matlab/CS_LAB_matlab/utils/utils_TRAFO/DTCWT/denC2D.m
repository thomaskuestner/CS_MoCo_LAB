function y = denC2D(x,T);

% % Example
% s1 = double(imread('st.tif'));
% s = s1(:,:,3);
% x = s + 20*randn(size(s));
% T = 40;
% y = denC2D(x,T);
% imagesc(y)
% colormap(gray)
% axis image
% sqrt(mean(mean((y-s).^2)))

[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
J = 4;
w = cplxdual2D(x,J,Faf,af);
I = sqrt(-1);
% loop thru scales:
for j = 1:J
    % loop thru subbands
    for s1 = 1:2
        for s2 = 1:3
            C = w{j}{1}{s1}{s2} + I*w{j}{2}{s1}{s2};
            C = soft(C,T);
            w{j}{1}{s1}{s2} = real(C);
            w{j}{2}{s1}{s2} = imag(C);
        end
    end
end
y = icplxdual2D(w,J,Fsf,sf);

