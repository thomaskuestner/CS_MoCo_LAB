function [mhi] = RandMask_rect(tolm, toln, m, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shiqian Ma
% Date : 09/05/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = zeros(m,n);


for i=1:m;
	for j=1:n;
		if ( ( abs(round(i-m/2)-(rand-0.5)*m/2)  < tolm ) &&  ( abs(round(j-n/2)-(rand-0.5)*n/2)  < toln ) )
			M(i,j) = 1;
		end
	end
end

cm = round(m/2);
cn = round(n/2);

% M(cm,:) =1;
% M(:,cn) =1;
 M(cm,cn) =1;


% upper half plane mask (not including origin)
Mh = zeros(m,n);
Mh = M;
% Mh(cm+2:m,:) = 0;
% Mh(cm+1,cn+1:n) = 0;


M = ifftshift(M);
mi = find(M);
Mh = ifftshift(Mh);
mhi = find(Mh);
