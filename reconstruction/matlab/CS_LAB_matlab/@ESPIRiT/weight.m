function wres = weight(MapOp,x)

wres = squeeze(MapOp.eigenVals(:,:,1,:)) .* x;
