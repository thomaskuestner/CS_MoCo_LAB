function res=jinc(x)
%res=jinc(x)

[tmp,err] = besselj(1,pi*x+eps);
res = tmp./(pi*x+eps);
if sum(err(:))
    disp('error');
end


