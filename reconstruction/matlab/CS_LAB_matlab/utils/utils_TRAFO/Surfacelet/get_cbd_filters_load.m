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
%%	get_cbd_filters_load.m
%%	
%%	First created: 08-22-05
%%	Last modified: 04-13-06
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function flt = get_cbd_filters_load(bo, mode)

exp = ['load cbd_filters_' num2str(bo) '_'];
if (mode == 'd')
    exp = [exp 'dec.mat'];
else
    exp = [exp 'rec.mat'];
end

eval(exp);


%% This is how we construct the filters
% flt = cell(1,4);
% [cfH0, cfH1, cfG0, cfG1] = polycbd('D0', bo);
% 
% if strcmp(mode, 'd')
%     f = cfH0.coeffs;
%     flt{1} = f(2:end, 2:end).';
%     flt{2} = fliplr(cfH0.center.' - [1 1]);
%     f = cfH1.coeffs;
%     flt{3} = f(2:end, 3:end).';
%     flt{4} = fliplr(cfH1.center.' - [1 2]);
% else
%     f = cfG0.coeffs;
%     flt{1} = f(3:end, 3:end).';
%     flt{2} = fliplr(cfG0.center.' - [2 2]);
%     f = cfG1.coeffs;
%     flt{3} = f(3:end, 2:end).';
%     flt{4} = fliplr(cfG1.center.' - [2 1]);
% end


%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.
    