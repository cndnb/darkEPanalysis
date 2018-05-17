%endCount = floor((stopFreq-startFreq)/(jump*(1/dataCut)))
function AMPOUT = fakeDarkEPanalysis(data, chunkSize, jump, startFreq, endCount,dataLength)
%Creates plotting array
ampFreq = ones(endCount,6);
%Creates error array
%ampVar = ones(endCount,7);
for count = 1:endCount

[BETA,COV] = specFreqPower(data,...
(startFreq+((count-1)*jump*(1/dataLength))),chunkSize); %Finds BETA for each frequency
ampFreq(count,:) = BETA;
%ampVar(count,:) = [(startFreq+((count-1)*jump*(1/rows(data)))),COV];

endfor
AMPOUT = ampFreq;
%AMPVAR = ampVar;
endfunction

%Error is not handled here, but in the standard deviation of the coherent
%average. This is why ampVar is commented

%!test
%! searchAmp = 1e-18;
%! count = 101;
%! startFreq = 1e-3;
%! stopFreq = 2e-2;
%! t1=1:10000; t1=t1';
%! sFreq = startFreq+((count-1)*(1/rows(t1)));
%! t2 = searchAmp.*sin((2*pi*sFreq).*t1);
%! fData = [t1,t2];
%! fEndC = floor((stopFreq-startFreq)/(1/rows(fData)));
%! testAmp = fakeDarkEPanalysis(fData,10,1,1e-3,fEndC);
%! tSpecAmp = testAmp(count,6); %This is the pure sine component at the same frequency as the input
%! assert (abs(tSpecAmp - searchAmp) < 2*eps)

%testAmp has the form [freq, B values from createSineComponents], where the pure sine component
%is the 5th term in createSineComponents => column 6 in testAmp.