%This function creates a weight value depending on the chi squared spread of the data
function weightVal = frequencyVariance(data,designX,inFreq,chunkLength,linearColumn)

if (nargin != 5)
  usage('weightVal = frequencyVariance(data,designX,inFreq,chunkLength,linearColumn)');
endif

%Resonant frequency of pendulum
f0 = 1.9338e-3;                                                                 
%Step size is chunk length periods of the frequency
stepSize = ceil(chunkLength*(1/inFreq));
%Number of whole chunks that fit evenly into the data
evenDivisible = floor(rows(data)/stepSize);
%Finds the remainder and subtracts from total to find end value for whole divisions
endValue = (evenDivisible*stepSize);
if (endValue == 0)
  stepSize = rows(data);
  endValue = rows(data);
endif

%Creates accumulation array for chi squared (weight) values
sineSTDev = ones(rows(data),1);

for counter = 0:stepSize:endValue
  %Prevents counter from stepping over the end of the array
  if (counter+stepSize > endValue)
    break;
  endif
    %Makes sure that constant and linear terms are orthogonal
    removeConstant = designX(counter+1:(counter+stepSize),:);
    if (linearColumn != 0)
      removeConstant(:,linearColumn) = removeConstant(:,linearColumn) .- (counter + 1);
    endif
    try   %Fit to find Chi^2
      [sChkBeta, sChkSigma, sChkR, sChkErr, sChkCov] = ols2(data(counter+1:(counter+stepSize),2),removeConstant);
      sineSTDev(counter+1:(counter+stepSize),1) = sChkSigma.*ones(stepSize,1);   
    catch %If X'*X becomes degenerate near resonance, removes those fit parameters
      inFreq
      fflush(stdout);
      onResonanceX = [removeConstant(:,1:6),removeConstant(:,9:columns(designX))];
      [sChkBeta, sChkSigma, sChkR, sChkErr, sChkCov] = ols2(data(counter+1:(counter+stepSize),2),onResonanceX);
      sineSTDev(counter+1:(counter+stepSize),1) = sChkSigma.*ones(stepSize,1);
    end_try_catch
  
endfor

if(endValue < rows(sineSTDev))
  %Makes the end points equal to the last analyzed point.
  sineSTDev(endValue:rows(sineSTDev),1) = sineSTDev(endValue).*ones(rows(sineSTDev)-endValue+1,1);
endif

%Returns the completed array of variances
weightVal = 1./(sineSTDev);


endfunction

%!test
%! b=1:10000;
%! b=b';
%! b=[b,randn(10000,1)]; %Data array time 1-10000, random y values
%! fitFreq = [1/1000, 1/100, 1/10];
%! for num = 1:rows(fitFreq) %Checks multiple frequncies
%!  varRes = frequencyVariance(b,createSineComponents(b(:,1),fitFreq(num)),fitFreq(num),50,0);
%!  for count = 1:rows(varRes)
%!    assert(varRes(count) != 0) %Checks that the weight is never zero
%!  endfor
%! endfor

%!test
%! b=1:10000;
%! b=b';
%! b=[b,randn(10000,1)]; %Data array time 1-10000, random y values
%! fitFreq = [1/1000, 1/100, 1/10];
%! for num = 1:rows(fitFreq)
%!  [B,olsSigma,R,Err,Cov] = ols2(b(:,2),createSineComponents(b(:,1),fitFreq(num))); %Finds chi square of fit
%!  varRes = frequencyVariance(b,createSineComponents(b(:,1),fitFreq(num)),fitFreq(num),rows(b)*fitFreq(num),0); %Chunk Length is equal to length of data
%!  fVSigma = varRes(1); %Takes first point, all points of fV should be the same
%!  assert (fVSigma == (1./olsSigma))
%! endfor




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OLD CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Specifies a round number based on input frequency for length of fit
%stepSize = ceil(chunkLength*(1/inFreq));
%modStep = stepSize;
%sineSTDev = ones(length(data(:,1)),1);
%for counter=1:stepSize:length(data(:,1))
%    %Makes sure step doesn't go over the end of the data
%    if (counter+stepSize > length(data(:,1)))
%      modStep = length(data(:,1))-counter;
%    endif
%    if (modStep >10)
%    %OLS fit to specific chunk
%    designX = genSineSeed(data(counter:(counter+modStep),1),inFreq);
%    [sChkBeta, sChkSigma, sChkR, sChkErr, sChkCov] = ols2(data(counter:(counter+modStep),2),designX);
%    %Prevents divide by zero errors
%    if (sChkSigma == 0)
%      sChkSigma = 1;
%    endif
%    %Assigns variance to all points in the chunk
%    sineSTDev(counter:(counter+modStep),1) = sChkSigma.*ones(modStep+1,1);
%    modStep=stepSize;
%    else
%    sineSTDev(counter:(counter+modStep),1) = 1.*ones(modStep+1,1);
%    modStep = stepSize;
%    endif
%    
%endfor