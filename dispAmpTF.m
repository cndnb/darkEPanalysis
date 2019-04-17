function [ampOut,errOut] = dispAmpTF(driftFix,frequencies,columnSelector,displayOut,seattleLat,seattleLong,compassDir,startTime)
%driftFix: cell array, each day of data is its own cell element in rows.
%frequencies: column vector of frequencies to be scanned.
%columnSelector is a 1x10 matrix, 0 in entry turns off column, 1 in entry uses column in fit.
%displayOut = 0 -> does not show output; = 1 -> displays output on command line.

	if (nargin != 8)
		error('[AMP,ERR] = dispAmpTF(driftFix,frequencies,columnSelector,displayOut,seattleLat,seattleLong,compassDir,startTime)');
	endif
	if (size(columnSelector) != [1,10])
		error('columnSelector must be 1x10 matrix');
	endif
  
  %Creates matrix to put beta values into correct column of valueStuff
  fitMat = diag(columnSelector);
  countDown = columns(fitMat);
  while(countDown > 0)
     if(!columnSelector(countDown))
     fitMat(:,countDown) = [];
   endif
   countDown = countDown - 1;
  endwhile
	
	f0 = 1.9338e-3;
	numBETAVal = sum(columnSelector);
	endCount = rows(frequencies);

	%Creates array to collect chunk values for mean/stdev
	valueStuff = zeros(endCount,10,rows(driftFix));
	compVar = zeros(endCount,numBETAVal,rows(driftFix));
	compOut = zeros(endCount,3,rows(driftFix));
	ampError = zeros(endCount,6,rows(driftFix));  

	startCount = 1;

	%Catches first frequency equal to zero
	if(frequencies(1) == 0 && columnSelector(2) == 1)
		for count = 1:rows(driftFix)
			if(columnSelector(1) == 0)
				valueStuff(1,1,count) = mean(driftFix{count,1}(:,2));
			else
				valueStuff(1,2,count) = mean(driftFix{count,1}(:,2));
			endif
		endfor
		startCount = 2;
  	endif
	%Catches last frequency equal to Nyquist
  	if(frequencies(end) == .5 && columnSelector(2) == 1)
		for count = 1:rows(driftFix)
			designX = cos(pi*driftFix{count,1}(:,1));
			[altB,altS,altERR,altR,altCOV] = ols2(driftFix{count,1}(:,2),designX);
			if(columnSelector(1) == 0)
				valueStuff(end,1,count) = altB;
			else
				valueStuff(end,2,count) = altB;
			endif
		endfor
		endCount = endCount - 1;
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
			
			designX = createSineComponents(driftFix{secCount,1}(:,1),frequencies(count),constantMultiples,columnSelector);
			lBTA = columns(designX);
			BETA = zeros(lBTA,1);
			COV = zeros(lBTA,lBTA);	
			try
			[BETA,SIGMA,R,ERR,COV] = ols2(driftFix{secCount,1}(:,2),designX);
			catch
				frequencies(count)
				fflush(stdout);
			end_try_catch

			%Adds data to each column in collection arrays
				valueStuff(count,:,secCount) = (fitMat*BETA)';
			compVar(count,1:columns(COV),secCount) = diag(COV)';
       		endfor
     	endfor
	
	%Performs weighted average using variances of OLS fit
	if (rows(driftFix) > 1)
		%Weighted average over different bin sizes
    		valAvg = mean(valueStuff,3);

    		%Takes stdev of central values to more accurately represent real error
		    ampError(:,1:columns(valueStuff)) = std(valueStuff,0,3)./size(valueStuff,3);
	else %If only one bin
    		valAvg = valueStuff;
    		ampError = zeros(endCount,columns(valueStuff));
  	endif
	
	%Adds weighted average values to be in complex amplitude format
  	compAvg = [valAvg(:,2) + i.*valAvg(:,1),valAvg(:,4) + i.*valAvg(:,3),valAvg(:,6) + i.*valAvg(:,5),valAvg(:,7),valAvg(:,8),valAvg(:,10) + i.*valAvg(:,9)];
  	errAvg  = [ampError(:,2) + i.*ampError(:,1),ampError(:,4) + i.*ampError(:,3),ampError(:,6) + i.*ampError(:,5),ampError(:,7),ampError(:,8),ampError(:,10) + i.*ampError(:,9)];
	
	%Returns
  	ampOut = compAvg;
	errOut = errAvg;

endfunction

%!test %Shows that reducing problem to sin cos fit matches fft. Somehow this broke
%! seattleLat = 0; seattleLong = 0; compassDir = 0; startTime = 0; 
%! t = (1:86400/100)';
%! columnSelector = [1,1,0,0,0,0,0,0,0,0];
%! fData = [t,stdnormal_rnd(rows(t),1)];
%! driftFix = cell(1,1); driftFix{1,1} = fData;
%! 
%! freqArray = (0:rows(t)/2)'./rows(t);
%! fftOut = (2/rows(t)).*fft(fData(:,2))(1:rows(t)/2 + 1,:);fftOut([1,end],:) = fftOut([1,end],:)./2;
%! 
%! olsOut = ones(rows(freqArray),1);
%! [ampOut,errOut] = dispAmpTF({fData},freqArray,columnSelector,0,seattleLat,seattleLong,compassDir,startTime);
%! olsOut = ampOut(:,1);
%!
%! oC = preCalcComponents(0,seattleLat,seattleLong,compassDir,startTime);
%! rO = abs(abs(olsOut./oC(1,1))./abs(fftOut) - 1);
%! assert(rO < 1e-9);
