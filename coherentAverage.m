function [ampFreq,ampError] = coherentAverage(t,data,dataDivisions,accArray)
  if (nargin != 3)
    usage("[avg,err] = coherentAverage(t,val,divisions,accArray)");
  endif
  
  betaValues = columns(data)-1;
  dataCut = floor((rows(data))/dataDivisions);
  valueStuff = ones(dataCut,columns(data)-1);
  ampFreq = accArray;
  ampError = ampFreq;
  i = 1;
  %Runs the fitter over each bin to find the amplitude at each frequency
  for j=0:(dataCut):((dataDivisions-1)*dataCut)
    %Sums each bin into one array
    sAmp = data(j+1:j+dataCut,:);
    valueStuff(:,betaValues*i-(betaValues-1):betaValues*i) = sAmp;
    i=i+1;
  endfor
  
  %Sums values over each bin and then averages for the mean
  for i=1:dataDivisions
    ampFreq(:,2:(1+betaValues)) = ampFreq(:,2:(1+betaValues)) + valueStuff(:,betaValues*i-(betaValues-1):betaValues*i);
  endfor
  ampFreq(:,2:(1+betaValues)) = ampFreq(:,2:(1+betaValues))./dataDivisions;
  
  %Sums (x-mean(x))^2 and then divides by N-1 takes the sqrt
  for j=1:dataDivisions
    ampError(:,2:(1+betaValues)) = ampError(:,2:(1+betaValues)) + ...
    (valueStuff(:,betaValues*i-(betaValues-1):betaValues*i).-ampFreq(:,2:(1+betaValues))).^2;
  endfor
  ampError(:,2:(1+betaValues)) = sqrt(ampError(:,2:(1+betaValues))./(dataDivisions-1));
endfunction

%!test
%! t=1:4; t=t';
%! num = [1,2,3,4]';
%! data = [t,num];
%! [avg,err] = coherentAverage(data,2,[1,0;2,0]);
%! assert(abs(avg(:,2) .- [2,3]') < eps)