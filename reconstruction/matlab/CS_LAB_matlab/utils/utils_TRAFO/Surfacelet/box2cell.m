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
%%	cell2box.m
%%	
%%	First created: 08-22-05
%%	Last modified: 04-13-06
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Yc = box2cell(Xb, level)

%   Convert the subband data stored in an N-D array to a cell array
%
%   Xb: an N-D array containing the subbands
%
%   level: decomposition level along each dimension.
%
%   Yc: a cell array containing all the subbands
%
%   See also cell2box.m


%% modify the level variable
level(level == -1) = 0;
level = 2 .^ level;

%% dimension of the problem
N = ndims(Xb);
szX = size(Xb);
number_of_sub = prod(level);
mult = szX ./ level;

%% allocate the output cell array
Yc = cell(number_of_sub, 1);
subs_cell = ind2subs_cell(level, mult, [1 : number_of_sub]');
subs_array = cell(N,1); 


for n = 1 : number_of_sub
    
    %% get the corresponding subs_array
    for k = 1 : N
        subs_array{k} = [1 : mult(k)] + subs_cell{k}(n);
    end
    
    Yc{n} = Xb(subs_array{:});
    
end

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.
