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

for counter = 0:stepSize:endValue
    if (counter+stepSize > endValue)
      break;
    endif
    designX = createSineComponents(data(counter+1:(counter+stepSize),1),inFreq);
    [sChkBeta, sChkSigma, sChkR, sChkErr, sChkCov] = ols2(data(counter+1:(counter+stepSize),2),designX);
    sineSTDev(counter+1:(counter+stepSize),1) = sChkSigma.*ones(stepSize,1);
endfor

%Makes the end points equal to the last analyzed point.
sineSTDev(endValue:rows(sineSTDev),1) = sineSTDev(endValue).*ones(rows(sineSTDev)-endValue+1,1);

%Returns the completed array of variances
weightVal = sineSTDev;


endfunction


%!test
%! b=[1,2,3,4,5,6,7,8,9,10];
%! b=b';
%! b=[b,rand(10,1)-(.5.*ones(10,1))];
%! b(:,1) = 10000.*b(:,1);
%![B,S,R,Err,Cov] = ols2(b(:,2),createSineComponents(b(:,1),pi));
%! varRes = frequencyVariance(b,pi,rows(b)*pi);
%! assert (varRes(1) == S);




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