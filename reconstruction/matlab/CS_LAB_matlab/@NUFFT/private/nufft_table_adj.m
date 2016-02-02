function Xk = nufft_table_adj(st, X)
%function Xk = nufft_table_adj(st, X)
% adjoint of table-based nufft interpolation.
% in
%	st		structure from nufft_init
%	X [M,?]		DTFT values
% out
%	Xk [*Kd,?]	DFT coefficients
% Copyright 2004-3-30, Jeff Fessler and Yingying Zhang, University of Michigan

dd = length(st.Kd);

% t = omega / gamma
tm = zeros(size(st.om));
for id=1:dd
	gam = 2*pi / st.Kd(id);
	tm(:,id) = st.om(:,id) / gam;
end

if size(X,1) ~= st.M
	error 'X size problem'
end

% adjoint of phase shift
if isvar('st.phase_shift') & ~isempty(st.phase_shift)
	X = X .* conj(st.phase_shift);
end

if isreal(X)
	X = complexify(X);
end

if dd == 1
	Xk = interp1_table_adj_mex(X, st.h{1}, ...
		int32(st.Jd), int32(st.Ld), tm, int32(st.Kd(1)));

elseif dd==2
	Xk = interp2_table_adj_mex(X, st.h{1}, st.h{2}, ...
		int32(st.Jd), int32(st.Ld), tm, int32(st.Kd));

elseif dd==3
	Xk = interp3_table_adj_mex(X, st.h{1}, st.h{2}, st.h{3}, ...
		int32(st.Jd), int32(st.Ld), tm, int32(st.Kd));

else
	error '> 3d not done'
end
