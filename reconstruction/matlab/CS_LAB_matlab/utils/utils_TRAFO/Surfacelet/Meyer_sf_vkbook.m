function theta = Meyer_sf_vkbook(x)

%   The smooth passband function for constructing Meyer filters
%   This one is from the VK book
%    
%   Yue Lu and Minh N. Do
%   First Created: 08-26-05
%   Last Revision: 08-26-05
%

theta = 3 * x .^ 2 - 2 * x .^ 3;
theta(x <= 0) = 0;
theta(x >= 1) = 1;

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.