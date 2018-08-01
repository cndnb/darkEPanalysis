function [FAMP,FERR,FPHASE] = ampToPower(compAvg,freqArray,kappa,f0,Q)  
	if(nargin != 5)
		usage('[FAMP,FERR] = ampToPower(compAvg,freqArray,kappa,f0,Q)');
	endif
  
	%Accumulation array
	ampMod = ones(rows(freqArray),3); %1-Z, 2-PerpGamma 3-ParaGamma
	errMod = ones(rows(freqArray),3);

	%Divides by transfer function to get power(frequency)
	ampMod = compAvg./transferFunction(freqArray,kappa,f0,Q);
  
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
%! testMatrix = [freq,3,3,4,4,3,4]; %Easy to calculated numbers in test matrix
%! tAmp = abs(5/transferFunction(freq,kappa,f0,Q)); %Expected (comparison) value of FAMP from test matrix
%! tErr = abs(sqrt(((3^2)/(3^2+4^2))*(3^2)+((4^2)/(3^2+4^2))*(4^2))/transferFunction(freq,kappa,f0,Q)); %Expected (comparison) value of FERR from test matrix
%! compAmp = [freq,tAmp,tAmp,tAmp];
%! compErr = [freq,tErr,tErr,tErr];
%! [outAmp,outErr] = ampToPower(testMatrix,testMatrix,kappa,f0,Q); %Actual function output of FAMP/FERR
%! assert(compAmp == outAmp(:,1:4))
%! %assert(compErr == outErr)

%!test
%! freq = 1:1000; freq = freq'; freq = 1./freq; %Generates column vector of frequencies
%! kappa = 2; f0 = 1.9e-3; Q = 500000;
%! testMatrix = ones(rows(freq),7);
%! testMatrix(:,1) = freq;
%! [outAmp,outErr] = ampToPower(testMatrix,testMatrix,kappa,f0,Q);
%! for count = 2:columns(outErr)
%!   assert(outAmp(:,count) = sqrt(2)./transferFunction(freq,kappa,f0,Q))
%!   %assert(outErr(:,count) = 1 ./transferFunction(freq,kappa,f0,Q))
%! endfor
