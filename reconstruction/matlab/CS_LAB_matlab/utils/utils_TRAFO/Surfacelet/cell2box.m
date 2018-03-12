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

function Yb = cell2box(Xc, level)

%   Convert the subband data stored in cell array to an N-D array
%
%   Xc: a cell array containing all the subbands
%
%   level: decomposition level along each dimension.
%
%   Yb: an N-D array containing all the subbands
%
%   See also box2cell.m
%


%% modify the level variable
level(level == -1) = 0;
level = 2 .^ level;

%% dimension of the problem
N = length(level);
mult = size(Xc{1});
szY = mult .* level;
number_of_sub = length(Xc);


%% allocate the output N-D array, which is either real or complex valued
%% according to the input Xc.
Yb = repmat(Xc{1}(1), szY);   
subs_cell = ind2subs_cell(level, mult, [1 : number_of_sub]');
subs_array = cell(N,1); 

for n = 1 : number_of_sub
    
    %% get the corresponding subs_array
    for k = 1 : N
        subs_array{k} = [1 : mult(k)] + subs_cell{k}(n);
    end
    
    Yb(subs_array{:}) = Xc{n};
    
end

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.
