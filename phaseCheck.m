
%freq = 9e-3;
%numDays = 10;
%reShape = ones(numDays,2);
%for count = 1:numDays
%	dayRange = ((count-1)*86400 + 1):count*86400;
%	[dZ,dPeX,dPaX] = createSineComponents(O(dayRange,1),freq);
%	[b1,s1,r1,err1,cov1] = ols2(O(dayRange,2),dZ);
%	reShape(count,:) = b1';
%endfor
%compOut = reShape(:,2) + i.*reShape(:,1);
%phase = atan2(reShape(:,2),reShape(:,1));
%compOut = compOut.*exp(i.*phase).*(-i);
%assert(angle(compOut),zeros(numDays,1),2*eps);

%numDays = 3;
%phaseZ = ones(rows(preAvgZ),numDays + 1);
%phaseZ(:,1) = freqArray;
%for count = 1:numDays
%	phaseZ(:,count + 1) = atan(preAvgZ(:,2,count)./preAvgZ(:,1,count));
%endfor
%
%for count = 1:numDays
%	figure(count + 4);
%	semilogx(phaseZ(:,1),phaseZ(:,count + 1));
%endfor
t = (1:500000)';
testData = [t,cos(2*pi*(9e-3).*t)];
if(!exist('d'))
	importfakeDarkEP
endif
driftFix = cell(2,2);
driftFix{1,1} = d(195000:235000,:);
dX1 = [ones(rows(driftFix{1,1}),1),(1:rows(driftFix{1,1}))',sin(omegaEarth.*driftFix{1,1}(:,1)),cos(omegaEarth.*driftFix{1,1}(:,1))];
[b1,s1,r1,err1,cov1] = ols2(driftFix{1,1}(:,2),dX1);
driftFix{1,1}(:,2) = driftFix{1,1}(:,2) - dX1*b1;
driftFix{2,1} = d(370000:410000,:);
dX2 = [ones(rows(driftFix{2,1}),1),(1:rows(driftFix{2,1}))',sin(omegaEarth.*driftFix{2,1}(:,1)),cos(omegaEarth.*driftFix{2,1}(:,1))];
[b2,s2,r2,err2,cov2] = ols2(driftFix{2,1}(:,2),dX2);
driftFix{2,1}(:,2) = driftFix{2,1}(:,2) - dX2*b2;
driftFix{1,2} = testData(195000:235000,:);
driftFix{2,2} = testData(370000:410000,:);

testLength = cell2mat(driftFix(:,1));
fullLength = rows(testLength);
stopFreq = 1e-2;
startFreq = 1e-3;
endCount = floor((stopFreq-startFreq)/(1/fullLength))+1;

freqArray = ones(endCount,1);
for count = 1:endCount
	freqArray(count,1) = (startFreq+((count-1)*(1/fullLength)));
endfor

[compAvg,compOut] = dispAmpTF(driftFix(:,1),freqArray,0,0,1,0);
[tcA,tcO] = dispAmpTF(driftFix(:,2),freqArray,0,0,1,0);

%assert(angle(tcO(:,1,1)),angle(compOut(:,1,1)))
%assert(angle(tcO(:,1,2)),angle(compOut(:,1,2)))
%compOut = compOut.*exp(-i.*angle(compOut));

figure(1);
semilogx(freqArray,[angle(compOut(:,1,1)),angle(compOut(:,1,2))])
figure(2);
semilogx(freqArray,[angle(tcO(:,1,1)),angle(tcO(:,1,2))])
%assert(angle(tcO(:,1,1)),angle(compOut(:,1,1)))
%assert(angle(tcO(:,1,2)),angle(compOut(:,1,2)))

%testAvg = ones(endCount,3);
%for freq = 1:endCount
%	allBETA = ones(1,3,rows(driftFix));
%	for count = 1:rows(driftFix)
%		designX = createSineComponents(driftFix{count,1}(:,1),freqArray(freq));
%		[BETA,SIGMA,R,ERR,COV] = ols2(driftFix{count,1}(:,2),designX);
%		allBETA(1,:,count) = [BETA(2,1)+i.*BETA(1,1),BETA(4,1)+i.*BETA(3,1),BETA(6,1)+i.*BETA(5,1)];
%		allBETA(1,:,count) = allBETA(1,:,count).*exp(-i.*angle(allBETA(1,:,count)));
%	endfor
%	testAvg(freq,:) = mean(allBETA,3);
%endfor

%assert(compAvg,testAvg);
%assert(angle(compAvg),angle(testAvg));
