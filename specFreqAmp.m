function [retB,retC] = specFreqAmp(data,designX,weightVal)
  if (nargin != 3)
    usage('[BETA,COV] = specFreqAmp(data,designX,weightVal)');
  endif
  %This gives the beta values of the fitted output to the search frequency with weights.
  try
    [freqBeta, freqCov] = weightedOLS(data(:,2),designX,weightVal);
  catch
    onResonanceX = [designX(:,1:6),designX(:,9:columns(designX))];
    [freqBeta, freqCov] = weightedOLS(data(:,2),onResonanceX,weightVal);
  end_try_catch

  %Returns amplitudes of cosine/sine components at searchFreq
  retB = freqBeta';
  retC = freqCov';
endfunction

%!test
%! t = 1:10000; t=t';
%! A = 1; f = randn;
%! data = [t, A.*sin((2*pi*f).*t)];
%! weightVal = ones(rows(t),1);
%! checkF = randn;
%! designX = [sin((2*pi*checkF).*t), cos((2*pi*checkF).*t)];
%! [oB,oS,oR,oERR,oCOV] = ols2(data(:,2),designX);
%! [retB,retC] = specFreqAmp(data,designX,weightVal);
%! assert(oB,retB');