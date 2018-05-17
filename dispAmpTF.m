function [AMP,ERR] = dispAmpTF(driftFix,stopFreq,startFreq,jump,fullLength,dataDivisions,chunkSize,numBETAVal,test)

  if (nargin != 9)
    usage('[AMP,ERR] = dispAmpTF(driftFix,stopFreq,startFreq,jump,fullLength,dataDivisions,chunkSize,numBETAVal,test) (test = 0 for normal operation, 1 for testing)');
  endif
  fullDataCut = fullLength
  dataCut = floor((rows(driftFix))/dataDivisions);

  %endCount is the #rows of frequency matrix
  %endCout = total frequency band divided by the smallest frequency jump
  %Integer so that it can be used for indexing
  endCount = floor((stopFreq-startFreq)/(jump*(1/fullDataCut)))+1;
  %Accumulation arrays
  %Amplitude at each frequency
  ampFreq = zeros(endCount,numBETAVal + 1);
  %Error of each amplitude value
  ampError = zeros(endCount,numBETAVal + 1);

  %Assigns frequency values for the first column of the frequency and error arrays
  for count = 1:endCount
    ampFreq(count,1) = (startFreq+((count-1)*jump*(1/fullDataCut)));
  endfor
  ampError(:,1) = ampFreq(:,1);

  %Creates array to collect chunk values for mean/stdev
  valueStuff = zeros(endCount,numBETAVal*dataDivisions);
  
  valCounter = 1;
  %Runs the fitter over each bin to find the amplitude at each frequency
  for secCount = 0:(dataCut):((dataDivisions-1)*dataCut)
    secCount
  
    sAmp = ones(endCount,numBETAVal);
  
    for count = 1:endCount
      if (test)
        sAmp(count,:) = ones(1,numBETAVal);
      else
      count
      fflush(stdout);
      [BETA,COV] = specFreqPower(driftFix(secCount+1:secCount+dataCut,:),ampFreq(count,1),chunkSize); %Finds BETA for each frequency
      sAmp(count,:) = BETA;
      endif
    endfor
  
    valueStuff(:,numBETAVal*valCounter-(numBETAVal - 1):numBETAVal*valCounter) = sAmp;
    valCounter = valCounter + 1;
  endfor

  %Sums values over each bin and then averages for the mean
  for count=1:dataDivisions
    ampFreq(:,2:(numBETAVal + 1)) = ampFreq(:,2:(numBETAVal+1)) + valueStuff(:,numBETAVal*count-(numBETAVal - 1):numBETAVal*count);
  endfor
  ampFreq(:,2:(numBETAVal + 1)) = ampFreq(:,2:(numBETAVal+1))./dataDivisions;

  %Sums (x-mean(x))^2 and then divides by N-1 takes the sqrt
  for count=1:dataDivisions
    ampError(:,2:(numBETAVal + 1)) = ampError(:,2:(numBETAVal + 1)) + (valueStuff(:,numBETAVal*count-(numBETAVal-1):numBETAVal*count).-ampFreq(:,2:(numBETAVal + 1))).^2;
  endfor
  ampError(:,2:(numBETAVal + 1)) = sqrt(ampError(:,2:(numBETAVal + 1))./(dataDivisions-1));
  
  %Returns
  AMP = ampFreq;
  ERR = ampError;
endfunction

%!test
%! t = 1:10000; t=t';
%! tAmp = 1e-15;
%! freq = pi;
%! fakeData = [t, tAmp.*cos((2*pi*freq).*t)];
%! fullDataCut = rows(fakeData);  
%! startFreq = 1e-4;
%! stopFreq = 1e-2; 
%! jump = 1; chunkSize = 50; dataDivisions = 5; numBETAVal = 1;
%! endCount = floor((stopFreq-startFreq)/(jump*(1/fullDataCut)))+1;
%! ampFreq = ones(endCount,1);
%! for i = 1:endCount
%!    ampFreq(i,1) = (startFreq+((i-1)*jump*(1/fullDataCut)));
%! endfor
%! [AMP,ERR] = dispAmpTF(fakeData,stopFreq,startFreq,jump,fullDataCut,dataDivisions,chunkSize, numBETAVal,1);
%! assert(ampFreq == AMP(:,1))
%! assert(AMP(:,2) == ones(rows(ampFreq),1))
%! assert(ERR(:,2) == 0)