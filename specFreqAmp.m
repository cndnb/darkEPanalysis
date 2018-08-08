function [retB,retC] = specFreqAmp(data,designX,weightVal)
  if (nargin != 3)
    usage('[BETA,COV] = specFreqAmp(data,designX,weightVal)');
  endif
  %This gives the beta values of the fitted output to the search frequency with weights.
  [freqBeta, freqCov] = weightedOLS(data(:,2),designX,weightVal);

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
