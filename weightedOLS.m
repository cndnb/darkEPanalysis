function [retB,retC] = weightedOLS(t,Y,sigma,f)
   if (nargin != 4)
    usage ("BETA = weightedOLS (t,Y,sigma,f)");
  endif
   X = createSineComponents(t,f);
  [nr, nc] = size (X);
  [ry, cy] = size (Y);
  if (nr != ry)
    error ("ols2: incorrect matrix dimensions");
  endif

  W = diag(1./(sigma(:,1).^2));
  Z = X' * W * X;

  BETA = inv (Z) * X' * W * Y;
  COV = inv(Z);
  %out = [t, X*BETA];
  retB = BETA;
  retC = COV;
endfunction