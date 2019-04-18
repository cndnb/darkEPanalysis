function rtn = equatorialTransferFunction(freqArray,omegaEarth,I,kappa,Q,isExternal)
	if(nargin != 6)
		error("rtn = equatorialTransferFunction(freqArray,omegaEarth,I,kappa,Q,isExternal)");
	endif
	
	wF = (1e-18)/(1.725e-14);
	omegaFreq = (2*pi).*freqArray;
	if(isExternal)
		denom = -I.*(omegaFreq .- omegaEarth).^2 .+ (i*(sqrt(kappa*I)/Q)).*(omegaFreq .- omegaEarth) .+ kappa.*ones(rows(omegaFreq),1);
		rtn = 1./(wF.*denom);
	else
		denom = -I.*(omegaFreq .- omegaEarth).^2 .+ (kappa*(1+(i/Q))).*ones(rows(omegaFreq),1);
		rtn = 1./(wF.*denom);
	endif
endfunction

%!test
%! kappa = 1; Q = 1; I = 1; omegaEarth = 1; freq = 1/(2*pi);
%! isExternal = 0;
%! assert(equatorialTransferFunction(freq,omegaEarth,I,kappa,Q,isExternal),.5*(1-i));
%! isExternal = 1;
%! assert(equatorialTransferFunction(freq,omegaEarth,I,kappa,Q,isExternal),1);

%!test
%! kappa = 0; Q = 1; I = 1; omegaEarth = 2*pi*(1/86164);
%! freq = rand(100,1);
%! eAmp = -1./((2*pi).*freq .- omegaEarth).^2;
%! for count = 0:1
%! 	assert(equatorialTransferFunction(freq,omegaEarth,I,kappa,Q,count),eAmp);
%! endfor
