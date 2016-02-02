% test_mecv2_reconstruction.m -
% image reconstruction using mirror-extended curvelets

if(1)
  N = 512;
  H = N/2;
  [a,b] = ndgrid(0:N-1);
  x = a+b-N;
  x = (x-min(min(x)))/max(max(x-min(min(x))));
  ns = ceil(log2(N) - 3);
  nag = 16; 
  ci = 128;
  
  c = mefcv2(x,N,N,ns,nag);
  
  cfs = [];
  for s=1:length(c);    for w=1:length(c{s});      cfs = [cfs; abs(c{s}{w}(:))];    end;  end
  cfs = sort(abs(cfs)); cfs = cfs(end:-1:1);
  val = cfs(ci);
  d = c;
  cnt = 0;
  for s=1:length(c)
    for w=1:length(c{s})
      cnt = cnt + sum(sum(abs(c{s}{w})>=val+eps));
      d{s}{w} = c{s}{w}.*(abs(c{s}{w})>=val+eps);
    end
  end
  fprintf(1, 'number of coefficients used = %d\n', cnt);
  y = meicv2(d,N,N,ns,nag);
  figure;  imagesc(real(y)); axis equal; axis tight; colormap gray; colorbar;
  title('Partial reconstruction with mirror-extended curvelets');
end


