%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	SurfBox-MATLAB (c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%	Yue Lu and Minh N. Do
%%
%%	Department of Electrical and Computer Engineering
%%	Coordinated Science Laboratory
%%	University of Illinois at Urbana-Champaign
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%	NXdiamondmapping.m
%%	
%%	First created: 04-20-05
%%	Last modified: 04-13-06
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = NXdiamondmapping(N, beta)

%% Get the mapping kernel for the fan-shaped frequency response
%
%   See also ......
%

if rem(N, 2) ~= 1
    error('N must be an odd integer so that we have a zero-phase mapping function.');
end

%% Use the kaiser window to truncate and smooth the ideal sinc sequences

%% w = window(@kaiser, N, beta);
w = calculateKaiserwindow(N, beta);

ind = [1 : N]' - (N+1)/2;
w = w .* sinc(ind / 2);
w(N+1) = 0; 

%% Get the mapping kernel
m = zeros(2*N-1);
[k1, k2] = meshgrid(-N+1:N-1, -N+1:N-1);
a = k1 + k2 + (N+1)/ 2;
a(a <= 0) = N + 1;
a(a > N) = N + 1;
a = w(a); 
b = k1 - k2 + (N+1)/ 2;
b(b <= 0) = N + 1;
b(b > N) = N + 1;
b = w(b);

m = a .* b;
%% If you want to check the frquency response of this mapping kernel.
% mapp = m;
% dispfreqd2(CRISPfilter(mapp, [1;1]));

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.
