function ret = resonanceVariance(data, chunkLength)
  %Resonant frequency of pendulum
  f0 = 1.9338e-3;                                                                 
  %Step size is chunk length periods of the frequency
  stepSize = ceil(chunkLength*(1/f0));
  %Number of whole chunks that fit evenly into the data
  evenDivisible = floor(rows(data)/stepSize);
  %Finds the remainder and subtracts from total to find end value for whole divisions
  endValue = (evenDivisible*stepSize);
  %If no chunks fit, use largest length
  if (endValue == 0)
    stepSize = rows(data);
    endValue = rows(data);
  endif
  
  %Initializes design matrix for resonant sine wave
  designX = [sin((2*pi*f0).*(1:stepSize)'), cos((2*pi*f0).*(1:stepSize)')];
  %Initializes accumulation array
  weightVal = ones(rows(data),1);
  %For each chunk calculate the chi squared and save it.
  for count = 1:evenDivisible
    [b,s,r,err,cov] = ols2(data(stepSize*(count-1) + 1:stepSize*count,2),designX);
    weightVal(stepSize*(count-1) + 1:stepSize*count,1) = s*ones(stepSize,1);
  endfor
  %If the analysis stops before the end of the array, extend last known value
  if (endValue < rows(data))
    weightVal(endValue:end) = weightVal(stepSize*evenDivisible,1)*ones(rows(data)-endValue,1);
  endif
  
  %Weights are equal to 1/chisquared
  ret = 1 ./ weightVal;
endfunction

%!test
%! b=1:10000;
%! b=b';
%! b=[b,randn(10000,1)]; %Data array time 1-10000, random y values
%! varRes = resonanceVariance(b,50);
%! for count = 1:rows(varRes)
%!  assert(varRes(count) != 0) %Checks that the weight is never zero
%! endfor
%! endfor

%!test
%! b=1:10000;
%! b=b';
%! b=[b,randn(10000,1)]; %Data array time 1-10000, random y values
%! f0 = 1.9338e-3;
%! [B,olsSigma,R,Err,Cov] = ols2(b(:,2),createSineComponents(b(:,1),fitFreq(num))); %Finds chi square of fit
%! varRes = resonanceVariance(b,rows(b)*fitFreq(num)); %Chunk Length is equal to length of data
%! fVSigma = varRes(1); %Takes first point, all points of fV should be the same
%! assert (fVSigma == (1./olsSigma))