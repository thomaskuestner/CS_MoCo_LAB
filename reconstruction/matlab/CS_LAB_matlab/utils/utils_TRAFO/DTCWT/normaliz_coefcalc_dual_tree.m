

N = 2^9; 
L = 6;
x = zeros(N);

%[Fsf,Faf] = FSfarras;
[Fsf,Faf]  = AntonB;
[sf,af] = dualfilt1;

W_zero = cplxdual2D(x, L, Faf, af);

figure(1), clf
nor = [];

for scale = 1:L
    no = length(W_zero{scale}{1}{1}{1})/2;
    for part = 1:2
	for dir = 1:2
		for dir1 = 1:3
    		W = W_zero; 
    		W{scale}{part}{dir}{dir1}(no,no) = 1;
		y = icplxdual2D(W, L, Fsf, sf);
		nor{scale}{part}{dir}{dir1} = sqrt(sum(sum(y.^2)));
		%sqrt(sum(sum(y.^2)))
		end
	end
    end
end

save nor_dualtree nor

