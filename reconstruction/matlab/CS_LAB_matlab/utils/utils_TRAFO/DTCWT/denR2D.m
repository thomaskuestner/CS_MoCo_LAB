function y = denR2D(x,T);

% % Example
% s1 = double(imread('st.tif'));
% s = s1(:,:,3);
% x = s + 20*randn(size(s));
% T = 10;
% y = denR2D(x,T);
% imagesc(y)
% colormap(gray)
% axis image
% sqrt(mean(mean((y-s).^2)))

[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
J = 4;
w = dualtree2D(x,J,Faf,af);
% loop thru scales:
for j = 1:J
    % loop thru subbands
    for s1 = 1:2
        for s2 = 1:3
            w{j}{s1}{s2} = soft(w{j}{s1}{s2},T);
        end
    end
end
y = idualtree2D(w,J,Fsf,sf);

