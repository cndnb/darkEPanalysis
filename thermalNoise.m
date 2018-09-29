function rtn = thermalNoise(freqArray,kappa,Q,T,f0,aCN,N,isExternal) %Autocollimator noise ~1e-9
	if (nargin != 8)
		error('out = thermalNoise(freqArray,kappa,Q,T,f0,aCN,N,isExternal)');
	endif

	%Constants/variable initialization
	kB = 1.38064852e-23;
	tN = ones(rows(freqArray),1);
	
	%Computes correct noise to type of damping
	if (isExternal)
		tN = sqrt(((aCN./transferFunction(freqArray,kappa,f0,Q,isExternal)).^2 + ((4*kB*T*kappa)./(Q*2*pi*f0)))/N);
	else
		tN = sqrt(((aCN./transferFunction(freqArray,kappa,f0,Q,isExternal)).^2 + ((4*kB*T*kappa)./((Q*2*pi).*freqArray)))/N);
	endif

	%Returns output
	rtn = tN;
endfunction
