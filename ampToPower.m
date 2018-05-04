function [FAMP,FERR] = ampToPower(ampFreq,ampErr,kappa,f0,Q)
  
  
  ampMod = ones(rows(ampFreq),5); %1-Time, 2-ParGamma 3-perpGamma 4-z 5-sum
  ampMod(:,1) = ampFreq(:,1); %Gets frequencies
  ampMod(:,2) = sqrt(ampFreq(:,[3,5])'*ampFreq(:,[3,5]));%ParGamma
  ampMod(:,3) = sqrt(ampFreq(:,[2,4])'*ampFreq(:,[2,4]));%PerpGamma
  ampMod(:,4) = sqrt(ampFreq(:,6:7)'*ampFreq(:,6:7)); %Z
  ampMod(:,5) = sqrt(ampFreq(:,2:7)'*ampFreq(:,2:7)); %Sum
  %We have sum in quadrature amplitudes for each direction
  
  FAMP = abs(ampMod(:,2:7)./transferFunction(ampMod(:,1),kappa,f0,Q));
  %Divides by transfer function to get power(frequency)
  
  errMod = ones(rows(ampErr),4);
  errMod(:,1) = ampErr(:,1);
  errMod(:,2) = sqrt(((ampFreq(:,3).^2)./(ampFreq(:,3).^2 .+ampFreq(:,5).^2)).*ampError(:,3).^2 ...
.+((ampFreq(:,5).^2)./(ampFreq(:,3).^2 .+ampFreq(:,5).^2)).*ampError(:,5).^2);
  errMod(:,3) = sqrt(((ampFreq(:,2).^2)./(ampFreq(:,2).^2 .+ampFreq(:,4).^2)).*ampError(:,2).^2 ...
.+((ampFreq(:,4).^2)./(ampFreq(:,2).^2 .+ampFreq(:,4).^2)).*ampError(:,4).^2);
  errMod(:,4) = sqrt(((ampFreq(:,6).^2)./(ampFreq(:,6).^2 .+ampFreq(:,7).^2)).*ampError(:,6).^2 ...
.+((ampFreq(:,7).^2)./(ampFreq(:,6).^2 .+ampFreq(:,7).^2)).*ampError(:,7).^2);
  %Sum error would be disgusting to calculated and not that useful
  
  FERR = abs(errMod(:,2:7)./transferFunction(errMod(:,1),kappa,f0,Q));
endfunction