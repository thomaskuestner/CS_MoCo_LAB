function x = partial_iFDCT(X, sY, sX, is_real, meta)
%Apply an inverse fast discrete curvelet transformation to X.
%The input image is of the size MxN
%The resulting pyramidial structure of the curvelet coefficients is reconstructed, using
%the meta data and filling neglected coefficients with 0. 

[M, N] = size(meta.stackedPos);
coef   = zeros(meta.wholeLength, 1);
if sX*sY < M*N
    %perform backward trafo into smaller image
    X = zpad(X, M, N);
    coef(meta.stackedPos) = X(:);
else
    coef(meta.stackedPos) = X(:);
end

%reconstruct pyramidial structure of the curvelet coefficients from stacked data:
recIndex = 1;
C 	 = meta.wholeTransform;
for s=1:length(C) 
  for w=1:length(C{s})
    [cM, cN] = size(C{s}{w}); %extract size of current coefficient matrix
    C{s}{w}  = reshape(coef(recIndex:recIndex+cM*cN-1), cM, cN); %reconstruct coefficient matrix
    recIndex = recIndex + cM*cN;
  end
end

x = ifdct_wrapping(C, is_real, M, N);
x = imresize(x, [sY, sX],'nearest');
end
