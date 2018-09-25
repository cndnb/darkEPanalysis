function ret = transferFunction(inFreq,kappa,f0,Q,isExternal)
	if (nargin != 5)
		error("tau(frequency) = transferFunction(frequency,kappa,resonanceFreq,Q,isExternal)");
	endif
	if (columns(inFreq)!= 1)
		error("transferFunction can only take a frequency column vector");
	endif

	%Initialization of array
	tF = ones(rows(inFreq),1);

	%Chooses external vs. internal damping function
	if (isExternal)
		tF = (1/kappa)./((1 .- (inFreq./f0).^2).+((sqrt(-1).*inFreq)./(Q*f0)));
	else
		tF = (1/kappa)./((1 .- (inFreq./f0).^2).+(sqrt(-1)./Q));
	endif	

	%Returns final output
	ret = tF;
endfunction

%!test
%! inFreq = ones(10,1);
%! assert(abs(transferFunction(inFreq,1,2,2,1) - ((6-2*sqrt(-1))/5).*ones(10,1))<eps)
%! assert(abs(transferFunction(inFreq,1,2,2,0) - (12/13-((8*sqrt(-1))/13)).*ones(10,1))<eps)

%!test
%! inFreq = [1/10,1/100,1/100]'; %Take some frequencies
%! kappa = pi; f0 = e; Q = sqrt(2);
%! fVal = [transferFunction(inFreq,kappa,f0,Q,1),transferFunction(inFreq,kappa,f0,Q,0)];
%! rVal = ones(rows(inFreq),2);
%! for ind = 1:rows(inFreq)
%!  rVal(ind,1) = (1/kappa)/((1-(inFreq(ind)/f0)^2)+((sqrt(-1)*inFreq(ind))/(Q*f0)));
%!  rVal(ind,2) = (1/kappa)/((1-(inFreq(ind)/f0)^2)+(sqrt(-1)/Q));
%! endfor
%! assert(fVal == rVal) %Makes sure each column is the right number
