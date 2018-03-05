function theta = rcos(x)

%   The raised cosine function
%
%   Yue Lu and Minh N. Do
%   First Created: 08-26-05
%   Last Revision: 08-26-05
%

theta = 0.5 * (1 - cos(pi*x));

theta(x <= 0) = 0;
theta(x >= 1) = 1;

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.