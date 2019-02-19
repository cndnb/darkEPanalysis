function [ampOut,errOut] = dispAmpTF(driftFix,frequencies,linearColumn,noRes,displayOut,seattleLat,seattleLong,compassDir,startTime)

	if (nargin != 9)
		error('[AMP,ERR] = dispAmpTF(driftFix,frequencies,linearColumn,noRes,displayOut,seattleLat,seattleLong,compassDir,startTime)');
	endif
  
	f0 = 1.9338e-3;
	numBETAVal = columns(createSineComponents(1,1,seattleLat,seattleLong,compassDir,startTime));
	endCount = rows(frequencies);

	%Creates array to collect chunk values for mean/stdev
	valueStuff = zeros(endCount,numBETAVal,rows(driftFix));
	compVar = zeros(endCount,numBETAVal,rows(driftFix));
	compOut = zeros(endCount,3,rows(driftFix));
	errOut = zeros(endCount,3,rows(driftFix));
	ampError = zeros(endCount,numBETAVal,rows(driftFix));  

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
	
	%Outputs total number of frequencies to scan	
	if (displayOut)
		endCount
	endif
	
	%Runs the fitter over each bin to find the amplitude at each frequency
	for secCount = 1:rows(driftFix)
		if (displayOut)
			secCount
		endif
		constantMultiples = preCalcComponents(driftFix{secCount,1}(:,1),seattleLat,seattleLong,compassDir,startTime);
		for count = startCount:endCount
			if (displayOut)
				count
				fflush(stdout);
			endif
			
			designX = createSineComponents(driftFix{secCount,1}(:,1),frequencies(count),constantMultiples);
:w
			if(noRes)
				designX = designX(:,1:numBETAVal - 2);
			endif
			if (abs(f0 - frequencies(count)) < 4*(frequencies(2,1) - frequencies(1,1)) && !noRes)
				designX = designX(:,1:numBETAVal - 2);
				if (!noRes)
					frequencies(count)
            				fflush(stdout);
          			endif
        		endif
        		if (linearColumn != 0)
          			%Prevents linear and constant term from becoming degenerate
          			designX(:,linearColumn) = designX(:,linearColumn) .- (driftFix{secCount,1}(1,1));
        		endif
			
			[BETA,SIGMA,R,ERR,COV] = ols2(driftFix{secCount,1}(:,2),designX);

			%Adds data to each column in collection arrays
			valueStuff(count,1:rows(BETA),secCount) = BETA';
			compVar(count,1:columns(COV),secCount) = diag(COV)';
       		endfor
     	endfor
	
	%Performs weighted average using variances of OLS fit
	if (rows(driftFix) > 1)
		%Weighted average over different bin sizes
    		valAvg = mean(valueStuff,3);

    		%Takes stdev of central values to more accurately represent real error
		ampError = std(valueStuff,0,3)./size(valueStuff,3);
	else %If only one bin
    		valAvg = valueStuff;
    		ampError = zeros(rows(compOut),columns(compOut));
  	endif
	
	%Adds weighted average values to be in complex amplitude format
  	compAvg = [valAvg(:,2) + i.*valAvg(:,1),valAvg(:,4) + i.*valAvg(:,3),valAvg(:,6) + i.*valAvg(:,5),valAvg(:,8) + i.*valAvg(:,7),valAvg(:,9),valAvg(:,10)];
	
	%Adds errors on real/imaginary components to find error of modulus
	modErr = sqrt((1 ./(valAvg(:,2).^2 + valAvg(:,1).^2)).*((valAvg(:,2).^2).*(ampError(:,2).^2)+(valAvg(:,1).^2).*(ampError(:,1).^2)));
 
	%Returns
  	ampOut = compAvg;
	errOut = modErr;

endfunction

%!test %Checks that each column is equal to specAmpFreq at that frequency
%! seattleLat = pi/4;
%! seattleLong = pi/4;
%! compassDir = pi/6;
%! startTime = 0;
%! t= 1:2*86164; t=t';
%! Amp = 1;
%! freq = (1/100);
%! chunkSize = 10;
%! fData = [t,Amp.*sin((2*pi*freq).*t)];
%! weightVal = resonanceVariance(fData,chunkSize);
%! fData = [fData,weightVal];
%! dataDivisions = cell(2,1);
%! dataDivisions{1,1} = fData(1:86164,:);
%! dataDivisions{2,1} = fData(86165:2*86164,:);
%! linearColumn = 0;
%! freqArray = (0:rows(fData)/2)'./rows(fData);
%! freqArray([1;end],:) = []; 
%!
%! [compAvg,compOut] = dispAmpTF(dataDivisions,freqArray,linearColumn,1,1,seattleLat,seattleLong,compassDir,startTime);
%! compareArray = zeros(rows(freqArray),6,rows(dataDivisions));
%! compareVar = zeros(rows(freqArray),6,rows(dataDivisions));
%! 	for secCount = 1:rows(dataDivisions)
%!  		for count = 1:rows(freqArray)
%!    			designX = createSineComponents(dataDivisions{secCount,1}(:,1),freqArray(count));
%!   			[BETA,SIGMA,R,ERR,COV] = ols2(dataDivisions{secCount,1}(:,2),designX);
%!			BETA = BETA';
%!			compareVar(count,:,secCount) = diag(COV)';
%!   			compareArray(count,:,secCount) = BETA;
%! 		endfor
%!	endfor
%!	fCA = sum(compareArray.*compareVar,3)./sum(compareVar,3);
%!	fCA = [fCA(:,2) + i.*fCA(:,1),fCA(:,4) + i.*fCA(:,3),fCA(:,6) + i.*fCA(:,5)];
%!	ccO = [compareArray(:,2,:) + i.*compareArray(:,1,:),compareArray(:,4,:) + i.*compareArray(:,3,:),compareArray(:,6,:) + i.*compareArray(:,5,:)];
%!	assert(fCA,compAvg);
%!	assert(ccO,compOut);
%! endfor

%!test
%! t = 1:86164; t=t';
%! Amp = 1e-16;
%! f = 2*pi*(9e-3);
%! fData = [t,Amp.*sin(f.*t)];
%! dX = [ones(rows(t),1),t];
%! [b,s,r,err,cov] = ols2(fData(:,2),dX);
%! fData(:,2) = fData(:,2) - dX*b;
%! dD = cell(1,1);
%! dD{1,1} = fData;
%!
%! startFreq = 1e-3;
%! stopFreq = 1e-2;
%! fullLength = rows(fData);
%! freqArray = 1:(floor(fullLength/2));
%! freqArray = [0,freqArray];
%! freqArray = freqArray';
%! freqArray = freqArray./fullLength;
%!
%! tempStart = freqArray - startFreq.*ones(rows(freqArray),1);
%! tempEnd = freqArray - stopFreq.*ones(rows(freqArray),1);
%! pastMinStart = Inf;
%! pastMinEnd = Inf;
%! minIndStart = 0;
%! minIndEnd = 0;
%! for count = 1:rows(freqArray)
%!  if (abs(tempStart(count)) < pastMinStart)
%!    pastMinStart = abs(tempStart(count));
%!    minIndStart = count;
%!  endif
%!  if (abs(tempEnd(count)) < pastMinEnd)
%!    pastMinEnd = abs(tempEnd(count));
%!    minIndEnd = count;
%!  endif
%! endfor
%! indStart = minIndStart;
%! indEnd = minIndEnd;
%! freqArray = freqArray(indStart:indEnd,1);
%!
%! [compAvg,errOut] = dispAmpTF(dD,freqArray,0,0,0);
%! fftOut = (2/rows(t)).*fft(fData(:,2));
%! fftOut = fftOut(1:(rows(fftOut)/2 + 1),:);
%! ratioPlot = abs(compAvg(:,1))./abs(fftOut(indStart:indEnd,:)) - 1;
%! assert(ratioPlot < 1e-5);

