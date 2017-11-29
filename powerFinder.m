function [retB,retC] = powerFinder(data,searchFreq,chunkSize,threshold)
  %errorFixData removes high uncertainty points and leaves those that are less noisy.
  errorFixData=removeHighError(data(:,1),data(:,2),...
  frequencyVariance(data,searchFreq,chunkSize),threshold);
  %This gives the beta values of the fitted output to the search frequency with weights.
  [freqBeta, freqCov] = weightedOLS(errorFixData(:,1),errorFixData(:,2),errorFixData(:,3),searchFreq);
  %data(:,1),data(:,2),...
  %frequencyVariance(data,searchFreq,chunkSize),searchFreq);
  
  %Extracts the diagonal of the covariance matrix to get a matrix of variances
  freqCov = diag(freqCov);
  %Returns amplitudes of cosine/sine components at searchFreq
  retB = freqBeta';
  retC = freqCov';
endfunction