function [X, meta] = partial_FDCT(x, is_real, finest, nbscales, nbangles_coarse)
%Apply a forward fast discrete curvelet transformation to x.
%The input image is of the size MxN
%The resulting pyramidial structure of the curvelet coefficients is sorted by size
%and truncated to the size MN. This way, the transformed image X and the original image x
%have same size. The original pyramidial structure will be reconstructed while the iFDCT is performed, filling
%neglected coefficients with 0.

[M, N] 	= size(x);
C 	= fdct_wrapping(x, is_real, finest, nbscales, nbangles_coarse);

stacked = [];
for s=1:length(C) 
  for w=1:length(C{s})
    stacked = [stacked; C{s}{w}(:)];
  end
end

[sorted, indicees] 	= sort(abs(stacked), 'descend');
wholeLength 		= M*N;
X                   = reshape(stacked(indicees(1:wholeLength)),M,N);
stackedPos          = reshape(indicees(1:wholeLength),M,N);

meta 	= struct('wholeTransform', {C}, 'wholeLength', wholeLength, 'stackedPos', stackedPos); % meta data needed while computing the iFDCT
end
