function e = cden2(s,x,t)

% % Example
% s1 = double(imread('st.tif'));
% s = s1(:,:,3);% x = s + 20*randn(size(s));
% t = 0:5:45;
% e = cden2(s,x,t);
% plot(t,e);

N = length(t);
for k = 1:N
    y = denC2D(x,t(k));
    e(k) = sqrt(mean(mean((y-s).^2)));
end

