function [rtn,compOut] = dispAmpTF(driftFix,frequencies,linearColumn,weighted,displayOut)

  if (nargin != 5)
    usage('[AMP,ERR] = dispAmpTF(driftFix,frequencies,linearColumn,fitIsWeighted,displayOut)');
  endif

  f0 = 1.9338e-3;
  numBETAVal = columns(createSineComponents(1,1));
  endCount = rows(frequencies);

  %Creates array to collect chunk values for mean/stdev
  valueStuff = zeros(endCount,12,rows(driftFix));
  compVar = zeros(endCount,12,rows(driftFix));
  compOut = zeros(endCount,3,rows(driftFix));
  errOut = zeros(endCount,3,rows(driftFix));
  startCount = 1;

  %Catches first frequency equal to zero
  if(frequencies(1) == 0)
	for count = 1:rows(driftFix)
		valueStuff(1,:,count) = [0,mean(driftFix{count,1}(:,2)),0,0,0,0];
	endfor
	startCount = 2;
  endif
  if(frequencies(end) == .5)
	for count = 1:rows(driftFix)
		designX = sin(pi*driftFix{count,1}(:,1));
		[altB,altS,altERR,altR,altCOV] = ols2(driftFix{count,1}(:,2),designX);
		valueStuff(end,:,count) = [altB,0,0,0,0,0];
	endfor
	endCount = endCount - 1;
  endif
		
  
  %Runs the fitter over each bin to find the amplitude at each frequency
    for secCount = 1:rows(driftFix)
      if (displayOut)
        secCount
      endif
      for count = startCount:endCount
        if (displayOut)
          count
          fflush(stdout);
        endif
        designX = createSineComponents(driftFix{secCount,1}(:,1),frequencies(count));
        if (linearColumn != 0)
          %Prevents linear and constant term from becoming degenerate
          designX(:,linearColumn) = designX(:,linearColumn) .- (driftFix{secCount,1}(1,1));
        endif
	if (abs(frequencies(count)-f0) < 3e-5)
	  designX = designX(:,1:10);
	endif
	allBETA = 0;
	allCov = 0;
	if (weighted)
		[BETA,COV] = specFreqAmp(driftFix{secCount,1}(:,1:2),designX,driftFix{secCount,1}(:,3));
		allBETA = BETA;
		allCov = diag(COV)';
	else
		[BETA,SIGMA,R,ERR,COV] = ols2(driftFix{secCount,1}(:,2),designX);
		allBETA = BETA';
		allCov = diag(COV)';
	endif

	%Adds data to each column in collection arrays
	valueStuff(count,1:columns(allBETA),secCount) = allBETA;
	compVar(count,1:columns(allCov),secCount) = allCov;
       endfor
     endfor
	
	%Weights smaller variance more heavily
	compVar = 1 ./ compVar;
	if (rows(driftFix) > 1) %This is only used in testing
		%Weighted average over different bin sizes
    		valAvg = sum(valueStuff.*compVar,3)./sum(compVar,3);
    		%Sums (x-mean(x))^2 and then divides by N-1 takes the sqrt
    		ampError = std(valueStuff,0,3); %0 makes std use denominator N-1
	else
    		valAvg = valueStuff;
    		ampError = zeros(rows(compOut),columns(compOut));
  	endif
  	compOut = [valueStuff(:,2,:)+ i.*valueStuff(:,1,:),valueStuff(:,4,:)+i.*valueStuff(:,3,:),valueStuff(:,6,:)+i.*valueStuff(:,5,:)];
  	compAvg = [valAvg(:,2) + i.*valAvg(:,1),valAvg(:,4) + i.*valAvg(:,3),valAvg(:,6) + i.*valAvg(:,5)];
  	
	%Returns
  	rtn = compAvg;

endfunction

%!test %Checks that each column is equal to specAmpFreq at that frequency
%! t= 1:10000; t=t';
%! Amp = 1;
%! freq = randn*(1/100);
%! chunkSize = 10;
%! fData = [t,Amp.*sin((2*pi*freq).*t)];
%! weightVal = resonanceVariance(fData,chunkSize);
%! fData = [fData,weightVal];
%! dataDivisions = cell(2,1);
%! dataDivisions{1,1} = fData(1:5000,:);
%! dataDivisions{2,1} = fData(5001:10000,:);
%! linearColumn = 0;
%! freqArray = (0:rows(fData)/2)'./rows(fData);
%! freqArray([1;end],:) = []; 
%!
%! for isWeighted = 0:1
%! 	[compAvg,compOut] = dispAmpTF(dataDivisions{1,1},freqArray,linearColumn,isWeighted,0);
%!
%! 	compareArray = ones(rows(freqArray),3);
%! 	for count = 1:rows(freqArray)
%!   		designX = createSineComponents(t,freqArray(count,1));
%!		BETA = ones(6,1);
%!   		if (isWeighted)
%!     			[BETA,COV] = specFreqAmp(fData(:,1:2),designX,fData(:,3));
%!   		else
%!     			[BETA,COV] = ols2(fData(:,2),designX);
%!			BETA = BETA';
%!   		endif
%! 		compareArray(count,:) = [BETA(1,2) + i.*BETA(1,1),BETA(1,4) + i.*BETA(1,3),BETA(1,6) + i.*BETA(1,5)];
%!		
%! 	endfor
%! 	assert(compAvg,compareArray);
%! endfor
%!
%! for isWeighted = 0:1
%! 	compAvg = dispAmpTF(dataDivisions,freqArray,linearColumn,isWeighted,0,0);
%!
%! 	compareArray = zeros(rows(freqArray),6,rows(dataDivisions));
%!      compareVar = zeros(rows(freqArray),6,rows(dataDivisions));
%! 	for secCount = 1:rows(dataDivisions)
%!  		for count = 1:endCount
%!    			designX = createSineComponents(dataDivisions{secCount,1}(:,1),freqArray(count));
%!    			if(isWeighted)
%!    				[BETA,COV] = specFreqAmp(dataDivisions{secCount,1}(:,1:2),designX,dataDivisions{secCount,1}(:,3));
%!    			else
%!    				[BETA,COV] = ols2(dataDivisions{secCount,1}(:,2),designX);
%! 				BETA = BETA';
%!    			endif
%!			compareVar(count,:,secCount) = diag(COV)';
%!    			compareArray(count,:,secCount) = BETA;
%!  		endfor
%! 	endfor
%! 	fCA = sum(compareArray.*compareVar,3)./sum(compareArray,3);
%! 	assert(compAvg,compareArray);
%! endfor

%!#
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

