function X = nufft_table_interp(st, Xk)
%function X = nufft_table_interp(st, Xk)
% table-based 1D and 2D nufft 
% in
%	st	structure	formed by nufft_init (through nufft_init_table)
%	Xk	[*Kd,L]		over-sampled DFT coefficients
% out
%	X	[M,L]		NUFFT values
% Copyright 2004-3-30, Jeff Fessler and Yingying Zhang, University of Michigan

dd = length(st.Kd);

% t = omega / gamma
tm = zeros(size(st.om));
for id=1:dd
	gam = 2*pi / st.Kd(id);
	tm(:,id) = st.om(:,id) / gam;
end

if size(Xk,1) ~= prod(st.Kd), error 'Xk size problem', end

% force Xk to be complex, since this is needed for pointers in the mex files.
if isreal(Xk), Xk = complex(Xk); end

if dd == 1
	X = interp1_table_mex(Xk, st.h{1}, ...
		int32(st.Jd), int32(st.Ld), tm);

elseif dd==2
	Xk = reshape(Xk, st.Kd);
	X = interp2_table_mex(Xk, st.h{1}, st.h{2}, ...
		int32(st.Jd), int32(st.Ld), tm);

elseif dd==3
	Xk = reshape(Xk, st.Kd);
	X = interp3_table_mex(Xk, st.h{1}, st.h{2}, st.h{3}, ...
		int32(st.Jd), int32(st.Ld), tm);

else
	error '> 3d not done'
end

% apply phase shift
if isvar('st.phase_shift') & ~isempty(st.phase_shift)
	X = X .* st.phase_shift;
end
