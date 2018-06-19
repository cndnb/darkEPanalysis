function [retB,retC] = specFreqAmp(data,designX,inFreq,chunkSize,linearColumn)
  if (nargin != 5)
    usage('[BETA,COV] = specFreqAmp(data,designX,searchFreq,chunkSize,linearColumn)');
  endif
  
  weightVal = frequencyVariance(data,designX,inFreq,chunkSize,linearColumn);
  %This gives the beta values of the fitted output to the search frequency with weights.
  try
    [freqBeta, freqCov] = weightedOLS(designX,data(:,2),weightVal);
  catch
    inFreq
    fflush(stdout);
    onResonanceX = [designX(:,1:6),designX(:,9:columns(designX))];
    [freqBeta, freqCov] = weightedOLS(onResonanceX,data(:,2),weightVal);
  end_try_catch

  %Returns amplitudes of cosine/sine components at searchFreq
  retB = freqBeta';
  retC = freqCov';
endfunction

%!test
%! t1=1:100000; 
%! t1=t1'; 
%! tAmp = 1e-17;
%! freq = [1/10, 1/100,1/1000];
%! t2 = zeros(length(t1),1);
%! for count = 1:length(freq)
%!  t2 = t2 + tAmp.*sin((2*pi*(freq(1,count))).*t1);
%! endfor
%! linearColumn = columns(createSineComponents(1,1)) - 1;
%! fData = [t1,t2];
%! for count = 1:length(freq)
%!  designX = createSineComponents(fData(:,1),freq(count));
%!  [B,C] = specFreqAmp(fData,designX,freq(1,count),50,linearColumn);
%!  specAmp = B(1,5); %This is the pure sine term in createSineComponents
%!  assert (abs(1-(tAmp / specAmp)) < .1)
%!  %assert (abs(B(1,6)) < 3*eps)
%! endfor