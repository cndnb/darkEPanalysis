function [FAMP,FERR,FPHASE] = ampToPower(compAvg,errAvg,freqArray,kappa,f0,Q,sampleInterval,torsionFiltered,isExternal)
	if(nargin != 9)
		error('[FAMP,FERR,FPHASE] = ampToPower(compAvg,errAvg,freqArray,kappa,f0,Q,sampleInterval,torsionFiltered)');
	endif
  
	%Accumulation array
	ampMod = ones(rows(freqArray),6); %1-Z, 2-PerpGamma 3-ParaGamma 4-drift 5-Constant offset 6-res freq
	errMod = ones(rows(freqArray),3);

	%Divides by transfer function to get power(frequency)
	if(torsionFiltered)
		ampMod = compAvg./transferFunction(freqArray,kappa,f0,Q,isExternal)./twoPointTransfer(freqArray,f0,sampleInterval);
	else
		ampMod = compAvg./transferFunction(freqArray,kappa,f0,Q,isExternal);
	endif
  
	%Return
	FAMP = [freqArray,abs(ampMod)];
  
	%Phase return
	FPHASE = [freqArray,angle(ampMod)];
	
	%Multiplies errors by same constant
	if(torsionFiltered)
		compErrAvg = errAvg./transferFunction(freqArray,kappa,f0,Q,isExternal)./twoPointTransfer(freqArray,f0,sampleInterval);
	else
		compErrAvg = errAvg./transferFunction(freqArray,kappa,f0,Q,isExternal);
	endif
	

	%Calculates error on modulus of calculated values, they should not be zero. If the real and imaginary components are both zero, they were artificially zeroed so the final error should also be artificially zeroed.
	try
		errMod(:,1) = sqrt((1 ./(real(compAvg(:,1)).^2 .+ imag(compAvg(:,1)).^2)).*((real(compAvg(:,1)).*real(compErrAvg(:,1))).^2 + (imag(compAvg(:,1)).*imag(compErrAvg(:,1))).^2));
	catch
		errMod(:,1) = zeros(rows(freqArray),1);
	end_try_catch
	try
		errMod(:,2) = sqrt((1 ./(real(compAvg(:,2)).^2 .+ imag(compAvg(:,2)).^2)).*((real(compAvg(:,2)).*real(compErrAvg(:,2))).^2 + (imag(compAvg(:,2)).*imag(compErrAvg(:,2))).^2));
	catch
		errMod(:,2) = zeros(rows(freqArray),1);
	end_try_catch
	try
		errMod(:,3) = sqrt((1 ./(real(compAvg(:,3)).^2 .+ imag(compAvg(:,3)).^2)).*((real(compAvg(:,3)).*real(compErrAvg(:,3))).^2 + (imag(compAvg(:,3)).*imag(compErrAvg(:,3))).^2));
	catch
		errMod(:,3) = zeros(rows(freqArray),1);
	end_try_catch
  
	errOut = abs(errMod);
	FERR = [freqArray,errOut];
endfunction

%!test
%! freqArray = [1/10,1/100,1/1000]; kappa = 2; f0 = 1.9e-3; Q = 500000; %Initialize random variables
%! sampleInterval = 1;
%! for count = 1:rows(freqArray)
%! freq = freqArray(count);
%! for isExternal = 0:1
%! 	for isTorsionFiltered = 0:1
%! 		testMatrix = [3+4i,3+4i,3+4i]; %Easy to calculate numbers in test matrix
%! 		tAmp = ones(1,3);
%! 		tErr = ones(1,3);
%! 		if(isTorsionFiltered)%Gets expected values of FAMP and FERR from test matrix
%!    			tAmp = abs((3+4i)/transferFunction(freq,kappa,f0,Q,isExternal)/twoPointTransfer(freq,f0,sampleInterval));
%!    			errVal = 1/transferFunction(freq,kappa,f0,Q,isExternal)/twoPointTransfer(freq,f0,sampleInterval);
%!    			tErr = sqrt((1/(25))*((3*real(errVal))^2 + (4*imag(errVal))^2));
%!   		else
%!    			tAmp = abs((3+4i)/transferFunction(freq,kappa,f0,Q,isExternal));
%!    			errVal = 1/transferFunction(freq,kappa,f0,Q,isExternal);
%!		    	tErr = sqrt((1/25)*((3*real(errVal))^2 + (4*imag(errVal))^2));
%!   		endif
%!   		compAmp = [freq,tAmp,tAmp,tAmp];
%!   		compErr = [freq,tErr,tErr,tErr];
%!   		[outAmp,outErr] = ampToPower(testMatrix,ones(1,3),freq,kappa,f0,Q,1,isTorsionFiltered,isExternal);
%!   		assert(compAmp,outAmp)
%!   		assert(compErr,outErr)
%!  	endfor
%! endfor
%! endfor
