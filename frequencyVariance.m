function ret = frequencyVariance(data,inFreq,chunkLength)

%Specifies a round number based on input frequency for length of fit
stepSize = ceil(chunkLength*(1/inFreq));
modStep = stepSize;
sineSTDev = ones(length(data(:,1)),1);
for counter=1:stepSize:length(data(:,1))
    %Makes sure step doesn't go over the end of the data
    if (counter+stepSize > length(data(:,1)))
      modStep = length(data(:,1))-counter;
    endif
    if (modStep >10)
    %OLS fit to specific chunk
    designX = genSineSeed(data(counter:(counter+modStep),1),inFreq);
    [sChkBeta, sChkSigma, sChkR, sChkErr, sChkCov] = ols2(data(counter:(counter+modStep),2),designX);
    %Prevents divide by zero errors
    if (sChkSigma == 0)
      sChkSigma = 1;
    endif
    %Assigns variance to all points in the chunk
    sineSTDev(counter:(counter+modStep),1) = sChkSigma.*ones(modStep+1,1);
    modStep=stepSize;
    else
    sineSTDev(counter:(counter+modStep),1) = 1.*ones(modStep+1,1);
    modStep = stepSize;
    endif
    
endfor

%Returns the completed array of variances
ret = sineSTDev;

endfunction