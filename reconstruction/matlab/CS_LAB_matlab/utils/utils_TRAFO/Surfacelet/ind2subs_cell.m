function subs_cell = ind2subs_cell(siz, mult, ndx)

%% This is a modification of the Matlab function ind2sub.m
%% Used in box2cell.m and cell2box.m
%%

siz = double(siz);

n = length(siz);
subs_cell = cell(n, 1);

k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
for i = n:-1:1,
  v = floor(ndx/k(i));  
  subs_cell{i} = v * mult(i); 
  ndx = rem(ndx,k(i));
end

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.   
