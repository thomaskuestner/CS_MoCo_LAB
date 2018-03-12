function [u0,mask]=cleanOF3D(u,mask)

% code modified/optimised:
% - mask in single precision
% 26.08.2015, Verena Neumann/Thomas Kuestner

% ::To do
% if nargin<2
%     mask=masking(u,'Laplacian',0.99);
% end

% generate the NaN mask
for i=1:3
    mask{i} = not(isnan(u{i}));
    mask{i} = single(mask{i});      %single/double
    badidx = find(mask{i}==0);
    u{i}(badidx)=0;
end

%%
% Diffusion
h(:,:,1) = [0, 0, 0; 0, 1, 0; 0, 0, 0];
h(:,:,2) = [0, 1, 0; 1, 0, 1; 0, 1, 0];
h(:,:,3) = [0, 0, 0; 0, 1, 0; 0, 0, 0];
h=convn(h,h);
u0=u;
encore=1;
mask0=mask;
tn = 0;
while encore
    tn = tn+1;
    %% ::to do, this number maybe not suitable
    %% intelligently filling
%     if tn>1000
%         break;
%     end
    
    un1=imfilter(mask0{1},h,'symmetric');
    f=imfilter(mask0{1}.*u0{1},h,'symmetric');
    n=find(mask{1}==0&un1~=0);
    u0{1}(n)=f(n)./un1(n);
    mask0{1}=double(un1~=0);
    
    un2=imfilter(mask0{2},h,'symmetric');
    f=imfilter(mask0{2}.*u0{2},h,'symmetric');
    n=find(mask{2}==0&un2~=0);
    u0{2}(n)=f(n)./un2(n);
    mask0{2}=double(un2~=0);
    
    un3=imfilter(mask0{3},h,'symmetric');
    f=imfilter(mask0{3}.*u0{3},h,'symmetric');
    n=find(mask{3}==0&un3~=0);
    u0{3}(n)=f(n)./un3(n);
    mask0{3}=double(un3~=0);
    
    encore=double(min(un1(:))==0);
end

% % Computation of averages
% window=2;
% u0=u;
% encore=1;
% todo=ones(1,length(badidx));
% while encore
%     window=2*window-1;
%     h=ones(window,window);
%     un=imfilter(mask,h,'symmetric');
%     f=imfilter(mask.*u,h,'symmetric');
%
%     todonow=todo.*(un(badidx)>=10)';
%     k0=find(todonow==1);
%     todo=todo.*(1-todonow);
%     encore=(sum(todo)>=1);
%
%     if ~isempty(k0)
%         u0(badidx(k0))=f(badidx(k0))./un(badidx(k0));
%     end
% end

% [M,N]=size(u);
% x=(1:M)'*ones(1,N);
% y=ones(M,1)*(1:N);
% good=find(mask==1);
% fx=TriScatteredInterp(x(good),y(good),real(u(good)));
% fy=TriScatteredInterp(x(good),y(good),imag(u(good)));
% u0(badidx)=fx(x(badidx),y(badidx))+i*fy(x(badidx),y(badidx));

%u0=zeroLap(u,mask);