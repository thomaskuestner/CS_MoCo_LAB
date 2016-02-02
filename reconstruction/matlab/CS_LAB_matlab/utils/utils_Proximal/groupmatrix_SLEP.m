function [ D, DT, ind, nodesize, waveS_l1, waveS_l12, waveS_SLEP] = groupmatrix_SLEP(n1,n2,wavelevel,waveletname_l1,waveletname_l12)
% creates the vectors and matrices needed for x=altra( v, n, ind, nodes) from SLEP package
% D: vector with indices according to a depth-first-search(tree)
% DT: inverse mapping of the wavelet coeffs
% ind: matrix with groupindices and weights
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

if wavelevel > 4
    wavelevel = 4; % max used wavelevel for SLEP
end;

image = zeros(n1,n2);
[~, waveS_l1] = wavedec2(image,wavelevel,waveletname_l1);
[~, waveS_l12] = wavedec2(image,wavelevel,waveletname_l12);
% make image bigger to cope uneven filter response
if ~strcmp(waveletname_l12,'db1')
    n1 = n1 + ceil(n1/3);
    n2 = n2 + ceil(n2/3);
end;
image = zeros(n1,n2);
[C, waveS_SLEP] = wavedec2(image,wavelevel,waveletname_l12);
n_coeff = size(C,2);

% single tree in required order:
cH = cell(1,4);
cH{1,1} = [85];
cH{1,2} = [63 84; 21 42];
cH{1,3} = [57 62 78 83; 47 52 68 73; 15 20 36 41; 5 10 26 31];
cH{1,4} = [55 56 60 61 76 77 81 82; 53 54 58 59 74 75 79 80; 45 46 50 51 66 67 71 72; 43 44 48 49 64 65 69 70; 13 14 18 19 34 35 39 40; 11 12 16 17 32 33 37 38; 3 4 8 9 24 25 29 30; 1 2 6 7 22 23 27 28];

%% adjust to regular tree:
flag_adjusted_n1 = false;
n1_SLEP = waveS_SLEP(1,2);
while ~flag_adjusted_n1
    if (n1_SLEP*2 <= waveS_SLEP(3,1))
        if (n1_SLEP*4 <= waveS_SLEP(4,1))
            if (n1_SLEP*8 <= waveS_SLEP(5,1))
                flag_adjusted_n1 = true;
                n1_SLEP = n1_SLEP + 1;
            end;
        end;
    end;
    n1_SLEP = n1_SLEP -1;
end;
flag_adjusted_n2 = false;
n2_SLEP = waveS_SLEP(1,2);
while ~flag_adjusted_n2
    if (n2_SLEP*2 <= waveS_SLEP(3,2))
        if (n2_SLEP*4 <= waveS_SLEP(4,2))
            if (n2_SLEP*8 <= waveS_SLEP(5,2))
                flag_adjusted_n2 = true;
                n2_SLEP = n2_SLEP + 1;
            end;
        end;
    end;
    n2_SLEP = n2_SLEP -1;
end;
% extend to whole wavelet space
for k = 1:wavelevel
    % extend in x dim.
    cH_start{1,k} = cH{1,k};
    % x_loop = floor(waveS_SLEP(k+1,2)/(2)^(k-1)-1);
    x_loop = n2_SLEP-1;
    for i = 1:x_loop
        cH{1,k} = [cH{1,k} cH_start{1,k}+85*i];
    end;
    % extend in y dim.
    cH_start2{1,k} = cH{1,k};
    % y_loop = floor(waveS_SLEP(k+1,1)/(2)^(k-1)-1);
    y_loop = n1_SLEP-1;
%     y_mod = mod(waveS_SLEP(k+1,1)/(2)^(k-1)-1,1);
    for v = 1:y_loop
        cH{1,k} = [cH_start2{1,k}+85*(x_loop+1)*v; cH{1,k}];
    end;
%     if (y_mod ~= 0)
%         cH{1,k} = [cH_start2{1,k}+85*(x_loop+1)*(y_loop+1); cH{1,k}];
%     end;
end;

%coeff_per_direction = (n_coeff - waveS_SLEP(1,1)*waveS_SLEP(1,2))/3;
coeff_per_direction = n1_SLEP * n2_SLEP * ( 1+2*2+4*4+8*8 );
for k = 1:wavelevel
    cH{1,k} = cH{1,k} + waveS_SLEP(1,1)*waveS_SLEP(1,2); % skip scale coeffs
    cV{1,k} = cH{1,k} + coeff_per_direction; % +nc-S(1,1)*S(1,2) 
    cD{1,k} = cH{1,k} + 2*coeff_per_direction; 
    
    if (size(cH{1,k},1) < waveS_SLEP(k+1,1))
        diff_1 = waveS_SLEP(k+1,1)- size(cH{1,k},1);
        cH{1,k}(size(cH{1,k},1)+1:size(cH{1,k},1)+diff_1,:) = 0;
        cV{1,k}(size(cH{1,k},1)+1:size(cH{1,k},1)+diff_1,:) = 0;
        cD{1,k}(size(cH{1,k},1)+1:size(cH{1,k},1)+diff_1,:) = 0;
    end;
    if (size(cH{1,k},2) < waveS_SLEP(k+1,2))
        diff_2 = waveS_SLEP(k+1,2)- size(cH{1,k},2);
        cH{1,k}(:,size(cH{1,k},2)+1:size(cH{1,k},2)+diff_2) = 0;
        cV{1,k}(:,size(cH{1,k},2)+1:size(cH{1,k},2)+diff_2) = 0;
        cD{1,k}(:,size(cH{1,k},2)+1:size(cH{1,k},2)+diff_2) = 0;
    end;
end;

D = zeros(1,n_coeff);
index_last = waveS_SLEP(1,1)*waveS_SLEP(1,2);
    D(1:index_last) = 1:index_last;
    %D(1:index_last) = zeros(1,index_last); % zeros for cA


% for i = 2:wavelevel+1
%     D(1+index_last:waveS_SLEP(i,1)*waveS_SLEP(i,2)+index_last) = reshape(cH{1,i-1},waveS_SLEP(i,1),waveS_SLEP(i,2));
%     D(1+index_last+waveS_SLEP(i,1)*waveS_SLEP(i,2):2*waveS_SLEP(i,1)*waveS_SLEP(i,2)+index_last) = reshape(cV{1,i-1},waveS_SLEP(i,1),waveS_SLEP(i,2));
%     D(1+index_last+2*waveS_SLEP(i,1)*waveS_SLEP(i,2):3*waveS_SLEP(i,1)*waveS_SLEP(i,2)+index_last) = reshape(cD{1,i-1},waveS_SLEP(i,1),waveS_SLEP(i,2));
%     index_last = index_last+3*waveS_SLEP(i,1)*waveS_SLEP(i,2);
% end;

%% readjust the following to replace the one above:
for i=2:wavelevel+1 % coeff.cH{2,wavelevel+2-i}(1:S(i,1),1:S(i,2)) indices to make sure that uneven coeffs have the right #
   D(1+index_last:waveS_SLEP(i,1)*waveS_SLEP(i,2)+index_last) = reshape(cH{1,i-1}(1:waveS_SLEP(i,1),1:waveS_SLEP(i,2)),1,waveS_SLEP(i,1)*waveS_SLEP(i,2));
   D(waveS_SLEP(i,1)*waveS_SLEP(i,2)+1+index_last:2*waveS_SLEP(i,1)*waveS_SLEP(i,2)+index_last) = reshape(cV{1,i-1}(1:waveS_SLEP(i,1),1:waveS_SLEP(i,2)),1,waveS_SLEP(i,1)*waveS_SLEP(i,2));
   D(2*waveS_SLEP(i,1)*waveS_SLEP(i,2)+1+index_last:3*waveS_SLEP(i,1)*waveS_SLEP(i,2)+index_last) = reshape(cD{1,i-1}(1:waveS_SLEP(i,1),1:waveS_SLEP(i,2)),1,waveS_SLEP(i,1)*waveS_SLEP(i,2));
   index_last = 3*waveS_SLEP(i,1)*waveS_SLEP(i,2)+index_last;
   % index_parent = 3*waveS_SLEP(i-1,1)*waveS_SLEP(i-1,2)+index_parent;
end;

DT = zeros(1,n_coeff);
end_regular_tree = 3 * coeff_per_direction + waveS_SLEP(1,1)*waveS_SLEP(1,2);
l = end_regular_tree;
for i = 1:n_coeff
    if D(1,i) == 0
        l = l+1;
        D(1,i) = l;
    end;
    DT(1,D(1,i)) = i;   
end;

% nodecount = (n_coeff - waveS_SLEP(1,1)*waveS_SLEP(1,2))/85;
nodecount = n1_SLEP * n2_SLEP * 3;

% test:
% E = [1:1:256*256];
% F = E(D);

%ind = ind(1,:) start, ind(2,:) end index, ind(3,:) weight
ind = zeros(3,85*nodecount+1);
ind_size = 85*nodecount;
ind(1,1:85) = [1 2 3 4 1 6 7 8 9 6 11 12 13 14 11 16 17 18 19 16 1 22 23 24 25 22 27 28 29 30 27 32 33 34 35 32 ...
    37 38 39 40 37 22 43 44 45 46 43 48 49 50 51 48 53 54 55 56 53 58 59 60 61 58 43 64 65 66 67 64 69 70 71 72 69 74 75 76 77 74 79 80 81 82 79 64 1];
ind(2,1:85) = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 ...
    37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85];

% no weight / invert weight / normal weight
ind(3,1:85) = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
%ind(3,1:85) = [1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1/sqrt(21) 1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1/sqrt(21) ...
%    1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1/sqrt(21) 1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1 1 1 1 1/sqrt(5) 1/sqrt(21) 1/sqrt(85)];
% ind(3,1:85) = [1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) sqrt(21) 1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) sqrt(21) ...
%     1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) sqrt(21) 1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) 1 1 1 1 sqrt(5) sqrt(21) sqrt(85)];

for v = 1:nodecount-1
    ind(1,1:(v+1)*85) = [ind(1,1:v*85)  v*85+ind(1,1:85)];
    ind(2,1:(v+1)*85) = [ind(2,1:v*85)  v*85+ind(2,1:85)];
    ind(3,1:(v+1)*85) = [ind(3,1:v*85)  ind(3,1:85)];
end;

ind(1,1:ind_size) = ind(1,1:ind_size) + waveS_SLEP(1,1)*waveS_SLEP(1,2);
ind(2,1:ind_size) = ind(2,1:ind_size) + waveS_SLEP(1,1)*waveS_SLEP(1,2);

% coeffs outside regular tree:
% vec_add = [end_regular_tree:n_coeff];

% root group:
ind(1,:) = [ind(1,1:ind_size) 1];
ind(2,:) = [ind(2,1:ind_size) n_coeff]; %ind_size+waveS_SLEP(1,1)*waveS_SLEP(1,2)
ind(3,:) = [ind(3,1:ind_size) 1]; %weight of root node : sqrt((ind_size+S(1,1)*S(1,2))) or 1 

nodesize = size(ind(1,:),2);
end