function [rtn,compOut,errOut] = dispAmpTF(driftFix,frequencies,linearColumn,weighted,displayOut,phaseFix)

  if (nargin != 6)
    usage('[AMP,ERR] = dispAmpTF(driftFix,frequencies,linearColumn,fitIsWeighted,displayOut,phaseFix)');
  endif

  numBETAVal = columns(createSineComponents(1,1));
  endCount = rows(frequencies);

  %Creates array to collect chunk values for mean/stdev
  compOut = zeros(endCount,3,rows(driftFix));
  errOut = zeros(endCount,3,rows(driftFix));
  startCount = 1;

  %Catches first frequency equal to zero
  if(frequencies(1) == 0)
	for count = 1:rows(driftFix)
		compOut(1,:,count) = [mean(driftFix{count,1}(:,2)),0,0];
	endfor
	startCount = 2;
  endif
  if(frequencies(end) == .5)
	for count = 1:rows(driftFix)
		designX = sin(pi*driftFix{count,1}(:,1));
		[altB,altS,altERR,altR,altCOV] = ols2(driftFix{count,1}(:,2),designX);
		compOut(end,:,count) = [i.*altB,0,0];
	endfor
	endCount = endCount - 1;
  endif
		
  
  %Runs the fitter over each bin to find the amplitude at each frequency
  if (weighted) %Performs weighted OLS fit
    for secCount = 1:rows(driftFix)
      if (displayOut)
        secCount
      endif
      for count = startCount:endCount
        if (displayOut)
          count
          fflush(stdout);
        endif
        [dZ,dPeX,dPaX] = createSineComponents(driftFix{secCount,1}(:,1),frequencies(count));
        %if (linearColumn != 0)
        %  %Prevents linear and constant term from becoming degenerate
        %  designX(:,linearColumn) = designX(:,linearColumn) .- (driftFix{secCount,1}(1,1));
        %endif
        data = driftFix{secCount,1};
        %Fits a data divison with the correct portion of the previously calculated design matrix
        [ZBETA,ZCOV] = specFreqAmp(data(:,1:2),dZ,data(:,3));
        compOut(count,1,secCount) = ZBETA(1,2) + i.*ZBETA(1,1);
        [PeXBETA,PeXCOV] = specFreqAmp(data(:,1:2),dPeX,data(:,3));
        compOut(count,2,secCount) = PeXBETA(1,2) + i.*PeXBETA(1,1);
        [PaXBETA,PaXCOV] = specFreqAmp(data(:,1:2),dPaX,data(:,3));
        compOut(count,3,secCount) = PaXBETA(1,2) + i.*PaXBETA(1,1);
      endfor
    endfor
  else %Performs unweighted OLS fit
    for secCount = 1:rows(driftFix)
      if (displayOut)
        secCount
      endif
  
      for count = startCount:endCount
        if (displayOut)
          count
          fflush(stdout);
        endif
        [dZ,dPeX,dPaX] = createSineComponents(driftFix{secCount,1}(:,1),frequencies(count));
        %designX = createSineComponents(driftFix{secCount,1}(:,1),frequencies(count));
        %if (linearColumn != 0)
        %  %Prevents linear and constant term from becoming degenerate
        %  designX(:,linearColumn) = designX(:,linearColumn) .- (driftFix{secCount,1}(1,1));
        %endif
        %Fits without weight the design matrix to the data
        %try
          [ZBETA,ZSIGMA,ZR,ZERR,ZCOV] = ols2(driftFix{secCount,1}(:,2),dZ);
          [PeXBETA,PeXSIGMA,PeXR,PeXERR,PeXCOV] = ols2(driftFix{secCount,1}(:,2),dPeX);
          [PaXBETA,PaXSIGMA,PaXR,PaXERR,PaXCOV] = ols2(driftFix{secCount,1}(:,2),dPaX);
	  %[BETA,SIGMA,R,ERR,COV] = ols2(driftFix{secCount,1}(:,2),designX);
          errOut(count,:,secCount) = [ZERR(2,1) + i.*ZERR(1,1),PeXERR(2,1) + i.* PeXERR(1,1),PaXERR(2,1) + i.*PaXERR(1,1)];
        %catch
          %frequencies(count)
	  %fflush(stdout);
          %noResonance = [designX(:,1:6),designX(:,9:numBETAVal)];
          %[BETA,SIGMA,R,ERR,COV] = ols2(driftFix{secCount,1}(:,2),...
          %noResonance);
		      %BETA = [BETA(1:6,:);zeros(2,columns(BETA));BETA(7:end,:)];
        %end_try_catch
	compOut(count,1,secCount) = ZBETA(2,1) + i.*ZBETA(1,1);
	compOut(count,2,secCount) = PeXBETA(2,1) + i.*PeXBETA(1,1);
	compOut(count,3,secCount) = PaXBETA(2,1) + i.*PaXBETA(1,1);
	%compOut(count,1,secCount) = BETA(2,1) + i.*BETA(1,1);
	%compOut(count,2,secCount) = BETA(4,1) + i.*BETA(3,1);
	%compOut(count,3,secCount) = BETA(6,1) + i.*BETA(5,1);
	if(phaseFix)
		for compCount = 1:3
        		compOut(count,compCount,secCount) = compOut(count,compCount,secCount)*exp(-i*angle(compOut(count,compCount,secCount)));
			assert(angle(compOut(count,compCount,secCount)),0,2*eps);
		endfor
	endif
      endfor
    endfor
  endif

  if (rows(driftFix) > 1) %This is only used in testing
    %Sums values over each bin and then averages for the mean
    compAvg = mean(compOut,3);
    %Sums (x-mean(x))^2 and then divides by N-1 takes the sqrt
    ampError = std(compOut,0,3); %0 makes std use denominator N-1
  else
    compAvg = compOut;
    ampError = zeros(rows(compOut),columns(compOut));
  endif
  
  %Returns
  rtn = compAvg;

endfunction

%!test %Checks that each column is equal to specAmpFreq at that frequency
%! for isWeighted = 0:1
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
%! endCount = rows(freqArray);
%!
%! compAvg = dispAmpTF(dataDivisions,freqArray,linearColumn,isWeighted,0,0); %isWeighted = 1; displayOutput = 0
%!
%! compareArray = ones(endCount,3);
%! for count = 1:endCount
%!   [dZ,dPeX,dPaX] = createSineComponents(t,freqArray(count,1));
%!   if (isWeighted)
%!     [BETA1,COV] = specFreqAmp(fData(:,1:2),dZ,fData(:,3));
%!     [BETA2,COV] = specFreqAmp(fData(:,1:2),dPeX,fData(:,3));
%!     [BETA3,COV] = specFreqAmp(fData(:,1:2),dPaX,fData(:,3));
%!   else
%!     [BETA1,COV] = ols2(fData(:,2),dZ);
%!     [BETA2,COV] = ols2(fData(:,2),dPeX);
%!     [BETA3,COV] = ols2(fData(:,2),dPaX);
%!     BETA1 = BETA1'; BETA2 = BETA2'; BETA3 = BETA3';
%!   endif
%!   compareArray(count,:) = [BETA1(1,2)+i.*BETA1(1,1),BETA2(1,2)+i.*BETA2(1,1),BETA3(1,2)+i.*BETA3(1,1)];
%! endfor
%! assert(compAvg,compareArray);
%! endfor

%!test %Checks that mean works
%! for isWeighted = 0:1
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
%! endCount = rows(freqArray);
%!
%! compAvg = dispAmpTF(dataDivisions,freqArray,linearColumn,isWeighted,0,0);
%!
%! compareArray = zeros(endCount,3);
%! for secCount = 1:rows(dataDivisions)
%!  for count = 1:endCount
%!    [dZ,dPeX,dPaX] = createSineComponents(dataDivisions{secCount,1}(:,1),freqArray(count));
%!    if(isWeighted)
%!    [BETA1,COV] = specFreqAmp(dataDivisions{secCount,1}(:,1:2),dZ,dataDivisions{secCount,1}(:,3));
%!    [BETA2,COV] = specFreqAmp(dataDivisions{secCount,1}(:,1:2),dPeX,dataDivisions{secCount,1}(:,3));
%!    [BETA3,COV] = specFreqAmp(dataDivisions{secCount,1}(:,1:2),dPaX,dataDivisions{secCount,1}(:,3));
%!    else
%!    [BETA1,COV] = ols2(dataDivisions{secCount,1}(:,2),dZ);
%!    [BETA2,COV] = ols2(dataDivisions{secCount,1}(:,2),dPeX);
%!    [BETA3,COV] = ols2(dataDivisions{secCount,1}(:,2),dPaX);
%!    BETA1 = BETA1'; BETA2 = BETA2'; BETA3 = BETA3';
%!    endif
%!    tempComp = [BETA1(1,2) + i.*BETA1(1,1),BETA2(1,2)+i.*BETA2(1,1),BETA3(1,2)+i.*BETA3(1,1)];
%!    compareArray(count,:) = compareArray(count,:) + tempComp;
%!  endfor
%! endfor
%! compareArray = compareArray ./ rows(dataDivisions);
%! assert(compAvg,compareArray);
%! endfor

%!test
%! t = 1:10000; t=t';
%! Amp = 1e-16;
%! f = 2*pi*(randn/10);
%! fData = [t,Amp.*sin(f.*t)];
%! dX = [ones(rows(t),1),t];
%! [b,s,r,err,cov] = ols2(fData(:,2),dX);
%! fData(:,2) = fData(:,2) - dX*b;
%! dD = cell(1,1);
%! dD{1,1} = fData;
%! freqArray = ((0:(rows(t)/2))')./(rows(t));
%! [compAvg,errOut] = dispAmpTF(dD,freqArray,0,0,0,0);
%! fftOut = (1/rows(t)).*fft(fData(:,2));
%! fftOut = fftOut(1:(rows(fftOut)/2 + 1),:);
%! 
%! fakeSignal =(rows(t)/2).*ifft(-[compAvg(:,1);flip(conj(compAvg(2:end-1,1)))]);
%!
%! assert(abs(real(fakeSignal) .- fData(:,2)) < 4*eps)


