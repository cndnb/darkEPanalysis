function [FAMP,FERR,FPHASE] = ampToPower(compAvg,ampError,freqArray,kappa,f0,Q,sampleInterval,torsionFiltered,isExternal)
	if(nargin != 9)
		error('[FAMP,FERR,FPHASE] = ampToPower(compAvg,ampError,freqArray,kappa,f0,Q,sampleInterval,torsionFiltered)');
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
	
	%Calculates error on modulus of calculated values, they should not be zero. If the real and imaginary components are both zero, they were artificially zeroed so the final error should also be artificially zeroed.
	try
		errMod(:,1) = sqrt((1 ./(real(compAvg(:,1)).^2 .+ imag(compAvg(:,1)).^2)).*((real(compAvg(:,1)).*ampError(:,1)).^2 + (imag(compAvg(:,1)).*ampError(:,2)).^2));
	catch
		errMod(:,1) = zeros(rows(freqArray),1);
	end_try_catch
	try
		errMod(:,2) = sqrt((1 ./(real(compAvg(:,2)).^2 .+ imag(compAvg(:,2)).^2)).*((real(compAvg(:,2)).*ampError(:,3)).^2 + (imag(compAvg(:,2)).*ampError(:,4)).^2));
	catch
		errMod(:,2) = zeros(rows(freqArray),1);
	end_try_catch
	try
		errMod(:,3) = sqrt((1 ./(real(compAvg(:,3)).^2 .+ imag(compAvg(:,3)).^2)).*((real(compAvg(:,3)).*ampError(:,6)).^2 + (imag(compAvg(:,3)).*ampError(:,6)).^2));
	catch
		errMod(:,3) = zeros(rows(freqArray),1);
	end_try_catch
  
	errOut = abs(errMod./transferFunction(errMod(:,1),kappa,f0,Q,isExternal));
	FERR = [freqArray,errOut];
endfunction

%!test
%! freq = pi; kappa = 2; f0 = 1.9e-3; Q = 500000; %Initialize random variables
%! isExternal = 0;
%! testMatrix = [3+4i,3+4i,3+4i]; %Easy to calculate numbers in test matrix
%! tAmp = abs((3+4i)/transferFunction(freq,kappa,f0,Q,isExternal));%Expected (comparison) value of FAMP from test matrix
%! tErr = abs(sqrt(2)/transferFunction(freq,kappa,f0,Q,isExternal));%Expected (comparison) value of FERR from test matrix
%! compAmp = [freq,tAmp,tAmp,tAmp];
%! compErr = [freq,tErr,tErr,tErr];
%! [outAmp,outErr] = ampToPower(testMatrix,ones(1,3),freq,kappa,f0,Q,1,0,isExternal);%Actual function output of FAMP/FERR
%! assert(compAmp,outAmp)
%! %assert(compErr == outErr) % FIX ERROR OUTPUT
%! isExternal = 1; %FIX THIS CASE
