%This function creates a weight value depending on the chi squared spread of the data
function weightVal = frequencyVariance(data,inFreq,chunkLength)
%Step size is chunk length periods of the frequency
stepSize = ceil(chunkLength*(1/inFreq));
%Number of steps that fit evenly into the chunk
evenDivisible = floor(rows(data)/stepSize);
%Finds the remainder and subtracts from total to find end value for whole divisions
endValue = rows(data) - (rows(data) - (evenDivisible*stepSize));

betaValues = ones(evenDivisible,3);
sineSTDev = ones(rows(data),1);
count = 0;
for counter = 1:stepSize:endValue
    if (counter+stepSize > endValue)
      break;
    endif
    count = count+1;
    designX = createSineComponents(data(counter:(counter+stepSize),1),inFreq);
    [sChkBeta, sChkSigma, sChkR, sChkErr, sChkCov] = ols2(data(counter:(counter+stepSize),2),designX);
    sineSTDev(counter:(counter+stepSize),1) = sChkSigma.*ones(stepSize+1,1);
endfor
sineSTDev(endValue:rows(sineSTDev),1) = sineSTDev(endValue).*ones(rows(sineSTDev)-endValue+1,1);

%Returns the completed array of variances
weightVal = sineSTDev;


endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OLD CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



