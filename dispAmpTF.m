function [AMP,ERR] = dispAmpTF(driftFix,frequencies,endCount,chunkSize,numBETAVal,linearColumn,weighted,displayOut)

  if (nargin != 8)
    usage('[AMP,ERR] = dispAmpTF(driftFix,frequencies,endCount,chunkSize,numBETAVal,linearColumn,fitIsWeighted,displayOut)');
  endif

  
  %Accumulation arrays
  %Amplitude at each frequency
  ampFreq = zeros(endCount,numBETAVal);
  %Error of each amplitude value
  ampError = zeros(endCount,numBETAVal);


  %Creates array to collect chunk values for mean/stdev
  valueStuff = zeros(endCount,numBETAVal,rows(driftFix));
  %Runs the fitter over each bin to find the amplitude at each frequency
  if (weighted) %Performs weighted OLS fit
    for secCount = 1:rows(driftFix)
      if (displayOut)
        secCount
      endif
      weightVal = resonanceVariance(driftFix{secCount,1},chunkSize);
      for count = 1:endCount
        if (displayOut)
          count
          fflush(stdout);
        endif
        designX = createSineComponents(driftFix{secCount,1}(:,1),frequencies(count));
        if (linearColumn != 0)
          %Prevents linear and constant term from becoming degenerate
          designX(:,linearColumn) = designX(:,linearColumn) .- (driftFix{secCount,1}(1,1));
        endif
        %Fits a data divison with the correct portion of the previously calculated design matrix
        [BETA,COV] = specFreqAmp(driftFix{secCount,1},...
        designX,weightVal);
        valueStuff(count,:,secCount) = BETA;
      endfor
    endfor
  else %Performs unweighted OLS fit
    for secCount = 1:rows(driftFix)
      if (displayOut)
        secCount
      endif
  
      for count = 1:endCount
        if (displayOut)
          count
          fflush(stdout);
        endif
        designX = createSineComponents(driftFix{secCount,1}(:,1),frequencies(count));
        if (linearColumn != 0)
          %Prevents linear and constant term from becoming degenerate
          designX(:,linearColumn) = designX(:,linearColumn) .- (driftFix{secCount,1}(1,1));
        endif
        %Fits without weight the design matrix to the data
        try
          [BETA,SIGMA,R,ERR,COV] = ols2(driftFix{secCount,1}(:,2),...
          designX);
        catch
          noResonance = [designX(:,1:6),designX(:,9:numBETAVal)];
          [BETA,SIGMA,R,ERR,COV] = ols2(driftFix{secCount,1}(:,2),...
          noResonance);
		BETA = [BETA(1:6,:);zeros(2,columns(BETA));BETA(7:end,:)];
        end_try_catch
	valueStuff(count,:,secCount) = BETA;
      endfor
    endfor
  endif
  
  if (rows(driftFix) > 1) %This is only used in testing
    %Sums values over each bin and then averages for the mean
    ampFreq = mean(valueStuff,3);


    %Sums (x-mean(x))^2 and then divides by N-1 takes the sqrt
    ampError = std(valueStuff,0,3); %0 makes std use denominator N-1
  else
    ampFreq = valueStuff;
    ampError = [];
  endif
  
  %Returns
  AMP = ampFreq;
  ERR = ampError;
endfunction

%!test %Checks that each column is equal to specAmpFreq at that frequency
%! t= 1:10000; t=t';
%! Amp = 1;
%! freq = randn*(1/100);
%! fData = [t,Amp.*sin((2*pi*freq).*t)];
%! startFreq = 1e-3;
%! stopFreq = 1e-2;
%! chunkSize = 50;
%! endCount = floor((stopFreq-startFreq)/(1/rows(t)))+1;
%! dataDivisions = cell(1,1);
%! dataDivisions{1,1} = fData;
%! numBETAVal = columns(createSineComponents(1,1));
%! linearColumn = 0;
%!
%! freqArray = ones(endCount,1);
%! for count = 1:endCount
%!   freqArray(count,1) = (startFreq+((count-1)*(1/rows(t))));
%! endfor
%!
%! [ampFreq,ampErr] = dispAmpTF(dataDivisions,freqArray,endCount,chunkSize,...
%! numBETAVal,linearColumn,1,0); %isWeighted = 1; displayOutput = 0
%!
%! compareArray = ones(endCount,numBETAVal);
%! weightVal = resonanceVariance(fData,chunkSize);
%! for count = 1:endCount
%!   [BETA,COV] = specFreqAmp(fData,createSineComponents(t,freqArray(count,1)),weightVal);
%!   compareArray(count,:) = BETA;
%! endfor
%! assert(ampFreq,compareArray);

%!test %Checks that mean works
%! t= 1:20000; t=t';
%! Amp = 1;
%! freq = randn*(1/100);
%! fData = [t,Amp.*sin((2*pi*freq).*t)];
%! startFreq = 1e-3;
%! stopFreq = 1e-2;
%! chunkSize = 50;
%! dataDivisions = cell(2,1);
%! dataDivisions{1,1} = fData(1:10000,:);
%! dataDivisions{2,1} = fData(10001:20000,:);
%! endCount = floor((stopFreq-startFreq)/(1/rows(t)))+1;
%! numBETAVal = columns(createSineComponents(1,1));
%! linearColumn = numBETAVal - 1;
%!
%! freqArray = ones(endCount,1);
%! for count = 1:endCount
%!   freqArray(count,1) = (startFreq+((count-1)*(1/rows(t))));
%! endfor
%!
%! [ampFreq,ampErr] = dispAmpTF(dataDivisions,freqArray,endCount,chunkSize,...
%! numBETAVal,linearColumn,1,0);%isWeighted = 1; displayOutput = 0
%!
%! compareArray = zeros(endCount,numBETAVal);
%! for secCount = 1:rows(dataDivisions)
%!  weightVal = resonanceVariance(dataDivisions{secCount,1},chunkSize);
%!  for count = 1:endCount
%!    removeConstant = createSineComponents(dataDivisions{secCount,1}(:,1),freqArray(count,1));
%!    removeConstant(:,linearColumn) = removeConstant(:,linearColumn) .- (dataDivisions{secCount,1}(1,1));
%!    [BETA,COV] = specFreqAmp(dataDivisions{secCount,1},...
%!    removeConstant,weightVal);
%!    compareArray(count,:) = compareArray(count,:) + BETA;
%!  endfor
%! endfor
%! compareArray = compareArray ./ rows(dataDivisions);
%! assert(ampFreq,compareArray);
