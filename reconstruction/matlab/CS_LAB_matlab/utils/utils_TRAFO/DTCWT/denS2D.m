function y = denS2D(x,T)

% % Example
% s1 = double(imread('st.tif'));
% s = s1(:,:,3);
% x = s + 20*randn(size(s));
% T = 35;
% y = denS2D(x,T);
% imagesc(y)
% colormap(gray)
% axis image

[af, sf] = farras;
J = 4;
w = dwt2D(x,J,af);
% loop thru scales:
for j = 1:J
    % loop thru subbands
    for s = 1:3
        w{j}{s} = soft(w{j}{s},T);
    end
end
y = idwt2D(w,J,sf);

