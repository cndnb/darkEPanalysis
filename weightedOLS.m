function [retB,retC] = weightedOLS(X,Y,sigma)
   if (nargin != 3)
    usage ("[BETA, COV] = weightedOLS (X,Y,sigma)");
  endif
  
  [nr, nc] = size (X);
  [ry, cy] = size (Y);
  if (nr != ry)
    error ("ols2: incorrect matrix dimensions");
  endif

  W = diag(1./(sigma(:,1).^2));
  Z = X' * W * X;

  BETA = inv (Z) * X' * W * Y;
  
  %This is the error on the amplitude fit parameters.
  COV = inv(Z);
  
  
  retB = BETA;
  retC = COV;
endfunction

%!test
%! testVal = [1,2,1;,4,2,1;2,3,Inf;4,3,Inf];
%! [freqBeta, freqCov] = weightedOLS(ones(rows(testVal),1),testVal(:,2),testVal(:,3));
%! assert (freqBeta == 2)