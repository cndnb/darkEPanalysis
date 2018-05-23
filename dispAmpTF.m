function [AMP,ERR] = dispAmpTF(driftFix,frequencies,endCount,dataDivisions,chunkSize,numBETAVal,test,weighted)

  if (nargin != 8)
    usage('[AMP,ERR] = dispAmpTF(driftFix,frequencies,endCount,dataDivisions,chunkSize,numBETAVal,test,fitIsWeighted) (test = 0 for normal operation, 1 for testing)');
  endif
  
  dataCut = floor((rows(driftFix))/dataDivisions);

  
  %Accumulation arrays
  %Amplitude at each frequency
  ampFreq = zeros(endCount,numBETAVal);
  %Error of each amplitude value
  ampError = zeros(endCount,numBETAVal);


  %Creates array to collect chunk values for mean/stdev
  valueStuff = zeros(endCount,numBETAVal,dataDivisions);
  
  %Runs the fitter over each bin to find the amplitude at each frequency
  if (weighted)
    %Performs weighted OLS fit
    for secCount = 0:(dataDivisions-1)
      secCount
  
      sAmp = ones(endCount,numBETAVal);
  
      for count = 1:endCount
        if (test)
          sAmp(count,:) = ones(1,numBETAVal);
        else
        count
        fflush(stdout);
        %Fits a data divison with the correct portion of the previously calculated design matrix
        [BETA,COV] = specFreqAmp(driftFix((secCount*dataCut)+1:(secCount*dataCut)+dataCut,:),...
        createSineComponents(driftFix((secCount*dataCut)+1:(secCount*dataCut)+dataCut,1),frequencies(count)),frequencies(count),chunkSize);
        sAmp(count,:) = BETA;
        endif
      endfor
  
      valueStuff(:,:,secCount + 1) = sAmp;
    endfor
  else
    %Performs unweighted OLS fit
    for secCount = 0:(dataDivisions-1)
      secCount
  
      sAmp = ones(endCount,numBETAVal);
  
      for count = 1:endCount
        if (test)
          sAmp(count,:) = ones(1,numBETAVal);
        else
        count
        fflush(stdout);
        [BETA,SIGMA,R,ERR,COV] = ols2(driftFix((secCount*dataCut)+1:(secCount*dataCut)+dataCut,2),...
        createSineComponents(driftFix((secCount*dataCut)+1:(secCount*dataCut)+dataCut,1),frequencies(count,1)));
        sAmp(count,:) = BETA;
        endif
      endfor
  
      valueStuff(:,:,secCount + 1) = sAmp;
    endfor
  endif
  
  %Sums values over each bin and then averages for the mean
  ampFreq = mean(valueStuff,3);


  %Sums (x-mean(x))^2 and then divides by N-1 takes the sqrt
  ampError = std(valueStuff,0,3); %0 makes std use denominator N-1
  
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
%! frequencies = abs(randn);
%! designX = createSineComponents(t,frequencies);
%! [AMP,ERR] = dispAmpTF(fakeData,designX,endCount,dataDivisions,chunkSize, numBETAVal,1,1);
%! assert(AMP == ones(endCount,1))
%! assert(ERR == zeros(endCount,1))