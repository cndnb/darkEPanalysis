%function [FAMP,FERR] = ampToPower(ampFreq,ampError,freqArray,kappa,f0,Q)
function [FAMP,FERR] = ampToPower(pAZ,pAPeX,pAPaX,freqArray,kappa,f0,Q)
  
  if(nargin != 7)
    usage('[FAMP,FERR] = ampToPower(pAZ,pAPeX,pAPaX,freqArray,kappa,f0,Q)');
  endif
  %if(columns(ampFreq) < 7)
  %  usage('ampFreq = [Frequency,ParGammaCos,ParGammaSine,PerpGammaCos,PerpGammaSine,Zcos,ZSine,...]');
  %endif
  %if(columns(ampError) < 7)
  %  usage('ampError = [Frequency,ParGammaCos,ParGammaSine,PerpGammaCos,PerpGammaSine,Zcos,ZSine,...]');
  %endif
  
  ampMod = ones(rows(freqArray),2); %1-Time, 2-ParGamma 3-perpGamma 4-z 5-sum
  ampMod(:,1) = freqArray; %Gets frequencies
  
  valDim = size(pAZ);
  preAvgZ = ones(rows(freqArray),valDim(3));
  for count = 1:valDim(3)
    preAvgZ = sqrt(pAZ(:,1,count).^2 + pAZ(:,2,count).^2);
    preAvgPerpX = sqrt(pAPeX(:,1,count).^2 + pAPeX(:,2,count).^2);
    preAvgParaX = sqrt(pAPaX(:,1,count).^2 + pAPaX(:,2,count).^2);
  endfor
  
  ampMod(:,2) = sum(preAvgZ,2)./columns(preAvgZ);
  ampMod(:,3) = sum(preAvgPerpX,2)./columns(preAvgPerpX);
  ampMod(:,4) = sum(preAvgParaX,2)./columns(preAvgParaX);
  %ampMod(:,2) = sqrt(ampFreq(:,1).^2 + ampFreq(:,2).^2);
  %ampMod(:,2) = sqrt(ampFreq(:,3).^2 + ampFreq(:,5).^2);%ParGamma
  %ampMod(:,3) = sqrt(ampFreq(:,2).^2 + ampFreq(:,4).^2);%PerpGamma
  %ampMod(:,4) = sqrt(ampFreq(:,6).^2 + ampFreq(:,7).^2); %Z
  %ampMod(:,5) = sqrt(ampFreq(:,2).^2 + ampFreq(:,3).^2 + ampFreq(:,4).^2 + ...
  %ampFreq(:,5).^2 + ampFreq(:,6).^2 + ampFreq(:,7).^2); %Sum
  %We have sum in quadrature amplitudes for each direction
  
  ampMod(:,2:4) = abs(ampMod(:,2:4)./transferFunction(ampMod(:,1),kappa,f0,Q));
  %ampMod(:,2:5) = abs(ampMod(:,2:5)./transferFunction(ampMod(:,1),kappa,f0,Q));
  FAMP = ampMod;
  %Divides by transfer function to get power(frequency)
  
  errMod = ones(rows(freqArray),4);
  %errMod(:,1) = ampError(:,1);
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
%! assert(compErr == outErr)

%!test
%! freq = 1:1000; freq = freq'; freq = 1./freq; %Generates column vector of frequencies
%! kappa = 2; f0 = 1.9e-3; Q = 500000;
%! testMatrix = ones(rows(freq),7);
%! testMatrix(:,1) = freq;
%! [outAmp,outErr] = ampToPower(testMatrix,testMatrix,kappa,f0,Q);
%! for count = 2:columns(outErr)
%!   assert(outAmp(:,count) = sqrt(2)./transferFunction(freq,kappa,f0,Q))
%!   assert(outErr(:,count) = 1 ./transferFunction(freq,kappa,f0,Q))
%! endfor
