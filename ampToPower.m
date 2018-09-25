function [FAMP,FERR,FPHASE] = ampToPower(compAvg,freqArray,kappa,f0,Q,sampleInterval,torsionFiltered,isExternal)  
	if(nargin != 8)
		error('[FAMP,FERR,FPHASE] = ampToPower(compAvg,freqArray,kappa,f0,Q,sampleInterval,torsionFiltered)');
	endif
  
	%Accumulation array
	ampMod = ones(rows(freqArray),6); %1-Z, 2-PerpGamma 3-ParaGamma 4-EarthModulaton 5-drift 6-Constant offset
	errMod = ones(rows(freqArray),6);

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

	%errMod(:,2) = sqrt((1 ./(ampFreq(:,3).^2 .+ ampFreq(:,5).^2)).*((ampFreq(:,3).*ampError(:,3)).^2 + (ampFreq(:,5).*ampError(:,5)).^2));
	%errMod(:,3) = sqrt((1 ./(ampFreq(:,2).^2 .+ ampFreq(:,4).^2)).*((ampFreq(:,2).*ampError(:,2)).^2 + (ampFreq(:,4).*ampError(:,4)).^2));
	%errMod(:,4) = sqrt((1 ./(ampFreq(:,6).^2 .+ ampFreq(:,7).^2)).*((ampFreq(:,6).*ampError(:,6)).^2 + (ampFreq(:,7).*ampError(:,7)).^2));
	%Sum error would be disgusting to calculate and not that useful
  
	%errMod(:,2:4) = abs(errMod(:,2:4)./transferFunction(errMod(:,1),kappa,f0,Q));
	FERR = errMod;
endfunction

%!test
%! freq = pi; kappa = 2; f0 = 1.9e-3; Q = 500000; %Initialize random variables
%! testMatrix = [3+4i,3+4i,3+4i]; %Easy to calculated numbers in test matrix
%! tAmp = abs((3+4i)/transferFunction(freq,kappa,f0,Q)); %Expected (comparison) value of FAMP from test matrix
%! tErr = abs(sqrt(((3^2)/(3^2+4^2))*(3^2)+((4^2)/(3^2+4^2))*(4^2))/transferFunction(freq,kappa,f0,Q)); %Expected (comparison) value of FERR from test matrix
%! compAmp = [freq,tAmp,tAmp,tAmp];
%! compErr = [freq,tErr,tErr,tErr];
%! [outAmp,outErr] = ampToPower(testMatrix,freq,kappa,f0,Q,1,0,1); %Actual function output of FAMP/FERR
%! assert(compAmp,outAmp)
%! %assert(compErr == outErr)
