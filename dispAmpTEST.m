
t= 1:10000; t=t';
 Amp = 1;
 freq = randn*(1/100);
 chunkSize = 10;
 fData = [t,Amp.*sin((2*pi*freq).*t)];
 weightVal = resonanceVariance(fData,chunkSize);
 fData = [fData,weightVal];
 dataDivisions = cell(2,1);
 dataDivisions{1,1} = fData(1:5000,:);
 dataDivisions{2,1} = fData(5001:10000,:);
 linearColumn = 0;
 freqArray = (0:rows(fData)/2)'./rows(fData);
 freqArray([1;end],:) = []; 

 isWeighted = 0
 	[compAvg,compOut] = dispAmpTF(dataDivisions,freqArray,linearColumn,isWeighted,0);

 	compareArray = zeros(rows(freqArray),6,rows(dataDivisions));
      	compareVar = zeros(rows(freqArray),6,rows(dataDivisions));
 	for secCount = 1:rows(dataDivisions)
  		for count = 1:rows(freqArray)
    			designX = createSineComponents(dataDivisions{secCount,1}(:,1),freqArray(count));
    			if(isWeighted)
    				[BETA,COV] = specFreqAmp(dataDivisions{secCount,1}(:,1:2),designX,dataDivisions{secCount,1}(:,3));
    			else
   				[BETA,COV] = ols2(dataDivisions{secCount,1}(:,2),designX);
				BETA = BETA';
   			endif
			compareVar(count,:,secCount) = diag(COV)';
   			compareArray(count,:,secCount) = BETA;
 		endfor
	endfor
	fCA = sum(compareArray.*compareVar,3)./sum(compareVar,3);
	fCA = [fCA(:,2) + i.*fCA(:,1),fCA(:,4) + i.*fCA(:,3),fCA(:,6) + i.*fCA(:,5)];
	ccO = [compareArray(:,2,:) + i.*compareArray(:,1,:),compareArray(:,4,:) + i.*compareArray(:,3,:),compareArray(:,6,:) + i.*compareArray(:,5,:)];
	assert(fCA,compAvg);
	assert(ccO,compOut);
 %endfor

