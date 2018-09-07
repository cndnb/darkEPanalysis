function ret = twoPointTransfer(inFreq,f0,Interval) %Update to work with discrete data
	if (nargin != 3)
		usage("tau(frequency) = twoPointTransfer(frequency,resonanceFreq,timingInterval)");
	endif
	if (columns(inFreq)!= 1)
		error("transferFunction can only take a frequency column vector");
	endif

	
	T0 = 1 / f0;
	fT0 = floor(T0/(4*Interval)); %Unitless floor function => floor(T0) = Interval*(fT0), ceil(T0) = Interval*(fT0 + 1)
	weight = (T0/(4*Interval)) - fT0; %1 - distance from actual value to floor, weight = 1 when T0 = fT0


	ret = (1-weight).*cos((2*pi)*Interval*fT0.*inFreq) + weight.*cos((2*pi)*Interval*(fT0 + 1).*inFreq);
endfunction

%!test
%! freq = ones(10,1);
%! T0 = 100;
%! f0 = 1/T0;
%! eVal = cos(pi/(2*f0));
%! cVal = twoPointTransfer(freq,f0,1);
%! assert(cVal,eVal.*ones(rows(freq),1));

%!test
%! T0 = 2.75;
%! f0 = 1/T0;
%! inFreq = randn/10;
%! cVal = twoPointTransfer(inFreq,f0,1);
%! eVal =.25*cos((pi/2)*floor(T0)*inFreq) + .75*cos((pi/2)*ceil(T0)*inFreq);
