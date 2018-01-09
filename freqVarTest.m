%How many periods of the specific frequency are included in error fit
chunkSize = 50;
%Multiples of smallest usable frequency between amplitude points
jump = 100;
%Start of frequency scan
startFreq = 1e-4;
%End frequency scan
stopFreq = 1e-1;
%How many points are included in the coherent average bins
dataCut = 100000;

%endCount is the #rows of frequency matrix
endCount = (stopFreq-startFreq)/(jump*(1/dataCut));
%Creates plotting array
ampFreq = zeros(endCount,7);
%Creates error array
ampError = zeros(endCount,7);

for i =1:endCount
  weights=frequencyVariance(driftFix,(startFreq+((i-1)*jump*(1/dataCut))),chunkSize);
  (startFreq+((i-1)*jump*(1/dataCut)))
  fflush(stdout);
  for j=1:length(driftFix)
    if (weights(j,1)==1)
      j
        fflush(stdout);
    endif
  endfor
endfor