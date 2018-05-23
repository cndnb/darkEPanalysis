function [retB,retC] = specFreqAmp(data,designX,inFreq,chunkSize)
  if (nargin != 4)
    usage('[BETA,COV] = specFreqAmp(data,designX,searchFreq,chunkSize)');
  endif
  
  weightVal = frequencyVariance(data,designX,inFreq,chunkSize);
  %This gives the beta values of the fitted output to the search frequency with weights.
  [freqBeta, freqCov] = weightedOLS(designX,data(:,2),weightVal);

  %Returns amplitudes of cosine/sine components at searchFreq
  retB = freqBeta';
  retC = freqCov';
endfunction

%!test
%! t1=1:10000; 
%! t1=t1'; 
%! tAmp = 1e-17;
%! freq = [1/10, 1/100,1/1000];
%! t2 = zeros(length(t1),1);
%! for count = 1:length(freq)
%! t2 = t2 + tAmp.*sin((2*pi*(freq(1,count))).*t1);
%! endfor
%! fData = [t1,t2];
%! for count = 1:length(freq)
%! [B,C] = specFreqPower(fData,freq(1,count),10);
%! specAmp = B(1,5); %This is the pure sine term in createSineComponents
%! tAmp/specAmp
%! assert (abs(1-(tAmp / specAmp)) < 2*eps)
%! %assert (abs(B(1,6)) < 2*eps)
%! endfor