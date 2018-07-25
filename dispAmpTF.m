%function [AMP,ERR] = dispAmpTF(driftFix,frequencies,endCount,linearColumn,weighted,displayOut)
function [preAvgZ,preAvgPerpX,preAvgParaX] = dispAmpTF(driftFix,frequencies,endCount,linearColumn,weighted,displayOut)


  if (nargin != 6)
    usage('[AMP,ERR] = dispAmpTF(driftFix,frequencies,endCount,linearColumn,fitIsWeighted,displayOut)');
  endif

  numBETAVal = columns(createSineComponents(1,1));
  
  %Accumulation arrays
  %Amplitude at each frequency
  ampFreq = zeros(endCount,numBETAVal);
  %Error of each amplitude value
  ampError = zeros(endCount,numBETAVal);


  %Creates array to collect chunk values for mean/stdev
  preAvgZ = zeros(endCount,2,rows(driftFix));
  preAvgPerpX = zeros(endCount,2,rows(driftFix));
  preAvgParaX = zeros(endCount,2,rows(driftFix));

  %Runs the fitter over each bin to find the amplitude at each frequency
  if (weighted) %Performs weighted OLS fit
    for secCount = 1:rows(driftFix)
      if (displayOut)
        secCount
      endif
      %weightVal = resonanceVariance(driftFix{secCount,1},chunkSize);
      for count = 1:endCount
        if (displayOut)
          count
          fflush(stdout);
        endif
        [dZ,dPeX,dPaX] = createSineComponents(driftFix{secCount,1}(:,1),frequencies(count));
        if (linearColumn != 0)
          %Prevents linear and constant term from becoming degenerate
          designX(:,linearColumn) = designX(:,linearColumn) .- (driftFix{secCount,1}(1,1));
        endif
        data = driftFix{secCount,1};
        %Fits a data divison with the correct portion of the previously calculated design matrix
        [ZBETA,ZCOV] = specFreqAmp(data(:,1:2),dZ,data(:,3));
        preAvgZ(count,:,secCount) = ZBETA;
        [PeXBETA,PeXCOV] = specFreqAmp(data(:,1:2),dPeX,data(:,3));
        preAvgPerpX(count,:,secCount) = PeXBETA;
        [PaXBETA,PaXCOV] = specFreqAmp(data(:,1:2),dPaX,data(:,3));
        preAvgParaX(count,:,secCount) = PaXBETA;
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
        [dZ,dPeX,dPaX] = createSineComponents(driftFix{secCount,1}(:,1),frequencies(count));
        if (linearColumn != 0)
          %Prevents linear and constant term from becoming degenerate
          designX(:,linearColumn) = designX(:,linearColumn) .- (driftFix{secCount,1}(1,1));
        endif
        %Fits without weight the design matrix to the data
        %try
          [ZBETA,ZSIGMA,ZR,ZERR,ZCOV] = ols2(driftFix{secCount,1}(:,2),...
          dZ);
          [PeXBETA,PeXSIGMA,PeXR,PeXERR,PeXCOV] = ols2(driftFix{secCount,1}(:,2)-dZ*ZBETA,...
          dPeX);
          [PaXBETA,PaXSIGMA,PaXR,PaXERR,PaXCOV] = ols2(driftFix{secCount,1}(:,2)-dZ*ZBETA,...
          dPaX);
        catch
          rank(designX)
          designX
          fflush(stdout);
          %noResonance = [designX(:,1:6),designX(:,9:numBETAVal)];
          %[BETA,SIGMA,R,ERR,COV] = ols2(driftFix{secCount,1}(:,2),...
          %noResonance);
		      %BETA = [BETA(1:6,:);zeros(2,columns(BETA));BETA(7:end,:)];
        end_try_catch
	      preAvgZ(count,:,secCount) = ZBETA';
        preAvgPerpX(count,:,secCount) = PeXBETA';
        preAvgParaX(count,:,secCount) = PaXBETA';
      endfor
    endfor
  endif

%  if (rows(driftFix) > 1) %This is only used in testing
%    %Sums values over each bin and then averages for the mean
%    ampFreq = mean(valueStuff,3);
%
%
%    %Sums (x-mean(x))^2 and then divides by N-1 takes the sqrt
%    ampError = std(valueStuff,0,3); %0 makes std use denominator N-1
%  else
%    ampFreq = valueStuff;
%    ampError = zeros(rows(valueStuff),columns(valueStuff));
%  endif
%  
%  %Returns
%  rtn = valueStuff;
%  %AMP = ampFreq;
%  %ERR = ampError;
endfunction

%!test %Checks that each column is equal to specAmpFreq at that frequency
%! t= 1:10000; t=t';
%! Amp = 1;
%! freq = randn*(1/100);
%! chunkSize = 50;
%! fData = [t,Amp.*sin((2*pi*freq).*t)];
%! weightVal = resonanceVariance(fData,chunkSize);
%! fData = [fData,weightVal];
%! startFreq = 1e-3;
%! stopFreq = 1e-2;
%! endCount = floor((stopFreq-startFreq)/(1/rows(t)))+1;
%! dataDivisions = cell(1,1);
%! dataDivisions{1,1} = fData;
%! linearColumn = 0;
%!
%! freqArray = ones(endCount,1);
%! for count = 1:endCount
%!   freqArray(count,1) = (startFreq+((count-1)*(1/rows(t))));
%! endfor
%!
%! [ampFreq,ampErr] = dispAmpTF(dataDivisions,freqArray,endCount,...
%! linearColumn,1,0); %isWeighted = 1; displayOutput = 0
%!
%! compareArray = ones(endCount,columns(createSineComponents(1,1)));
%! for count = 1:endCount
%!   [BETA,COV] = specFreqAmp(fData(:,1:2),createSineComponents(t,freqArray(count,1)),fData(:,3));
%!   compareArray(count,:) = BETA;
%! endfor
%! assert(ampFreq,compareArray);

%!test %Checks that mean works
%! t= 1:20000; t=t';
%! Amp = 1;
%! freq = randn*(1/100);
%! chunkSize = 50;
%! fData = [t,Amp.*sin((2*pi*freq).*t)];
%! weightVal = resonanceVariance(fData,chunkSize);
%! fData = [fData,weightVal];
%! startFreq = 1e-3;
%! stopFreq = 1e-2;
%! dataDivisions = cell(2,1);
%! dataDivisions{1,1} = fData(1:10000,:);
%! dataDivisions{2,1} = fData(10001:20000,:);
%! endCount = floor((stopFreq-startFreq)/(1/rows(t)))+1;
%! linearColumn = columns(createSineComponents(1,1)) - 1;
%!
%! freqArray = ones(endCount,1);
%! for count = 1:endCount
%!   freqArray(count,1) = (startFreq+((count-1)*(1/rows(t))));
%! endfor
%!
%! [ampFreq,ampErr] = dispAmpTF(dataDivisions,freqArray,endCount,...
%! linearColumn,1,0);%isWeighted = 1; displayOutput = 0
%!
%! compareArray = zeros(endCount,columns(createSineComponents(1,1)));
%! for secCount = 1:rows(dataDivisions)
%!  for count = 1:endCount
%!    removeConstant = createSineComponents(dataDivisions{secCount,1}(:,1),freqArray(count,1));
%!    removeConstant(:,linearColumn) = removeConstant(:,linearColumn) .- (dataDivisions{secCount,1}(1,1));
%!    [BETA,COV] = specFreqAmp(dataDivisions{secCount,1}(:,1:2),...
%!    removeConstant,dataDivisions{secCount,1}(:,3));
%!    compareArray(count,:) = compareArray(count,:) + BETA;
%!  endfor
%! endfor
%! compareArray = compareArray ./ rows(dataDivisions);
%! assert(ampFreq,compareArray);
