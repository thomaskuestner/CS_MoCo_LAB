function res = mtimes(a,b)


if a.adjoint
	res = adjD(b);

else
	res = D(b);

end




    
