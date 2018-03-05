function L=GetLevel(i,Ms,n) 
%return the level of i-th coefficient
col = floor(i/n);
row = rem(i,n);

dlevels = length(Ms);

p=n/(2^(dlevels-1)); 

for j=1:dlevels
    if (col<=p) && (row<=p)
        L=j;
        break;
    else
        p=2*p;
    end
end