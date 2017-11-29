function [retB,retC] = powerFinder(data,searchFreq,chunkSize)
  %This gives the beta values of the fitted output to the search frequency with weights.
  [freqBeta, freqCov] = weightedOLS(data(:,1),data(:,2),...
  frequencyVariance(data,searchFreq,chunkSize),searchFreq);
  freqCov = diag(freqCov);
  %Returns amplitudes of cosine/sine components at searchFreq
  retB = freqBeta';
  retC = freqCov';
endfunction