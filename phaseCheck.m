testLength = cell2mat(driftFix(:,1));
fullLength = rows(driftFix);
stopFreq = 1e-2;
startFreq = 1e-3;
endCount = floor((stopFreq-startFreq)/(1/rows(t)))+1;

freqArray = ones(endCount,1);
for count = 1:endCount
	freqArray(count,1) = (startFreq+((count-1)*(1/rows(t))));
endfor

compAvg = dispAmpTF(driftFix,freqArray,0,0,1);

testAvg = ones(endCount,3);
for freq = 1:endCount
	allBETA = ones(1,3,rows(driftFix));
	for count = 1:rows(driftFix)
		designX = createSineComponents(driftFix{count,1}(:,1),freqArray(freq));
		[BETA,SIGMA,R,ERR,COV] = ols2(driftFix{count,1}(:,2),designX);
		allBETA(1,:,count) = [BETA(2,1)+i.*BETA(1,1),BETA(4,1)+i.*BETA(3,1),BETA(6,1)+i.*BETA(5,1)];
		allBETA(1,:,count) = allBETA(1,:,count).*exp(-i.*angle(allBETA(1,:,count)));
	endfor
	testAvg(freq,:) = mean(allBETA,3);
endfor

assert(compAvg,testAvg);
assert(angle(compAvg),angle(testAvg));
