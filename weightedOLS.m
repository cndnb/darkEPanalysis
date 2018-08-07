function [retB,retC] = weightedOLS(Y,X,sigma)
	if (nargin != 3)
    		error("[BETA, COV] = weightedOLS (X,Y,sigma); %Sigma is a column vector");
  	endif
  
	[nr, nc] = size (X);
	[ry, cy] = size (Y);
  	if (nr != ry)
    		error ("ols2: incorrect matrix dimensions");
  	endif

  	W = diag(1 ./(sigma(:,1).^2));
  	Z = X' * W * X;

  	BETA = inv (Z) * X' * W * Y;
  
  	%This is the error on the amplitude fit parameters.
 	COV = inv(Z);
  
  	%Returns
  	retB = BETA;
  	retC = COV;
endfunction

%!test
%! testVal = [1,2,1;,4,2,1;2,3,Inf;4,3,Inf];
%! [freqBeta, freqCov] = weightedOLS(testVal(:,2),ones(rows(testVal),1),testVal(:,3));
%! assert (freqBeta == 2)

%!test
%! testVal = [1,2,1;4,2,1;2,3,1;4,3,1];
%! designMatrix = ones(rows(testVal),1);
%! [freqBeta, freqCov] = weightedOLS(testVal(:,2),designMatrix,testVal(:,3));
%! [oB,oS,oR,oERR,oCOV] = ols2(testVal(:,2),designMatrix);
%! assert (oB,freqBeta');
