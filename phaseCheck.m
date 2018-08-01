
freq = 9e-3;
numDays = 10;
reShape = ones(numDays,2);
for count = 1:numDays
	dayRange = ((count-1)*86400 + 1):count*86400;
	[dZ,dPeX,dPaX] = createSineComponents(O(dayRange,1),freq);
	[b1,s1,r1,err1,cov1] = ols2(O(dayRange,2),dZ);
	reShape(count,:) = b1';
endfor
compOut = reShape(:,2) + i.*reShape(:,1);
phase = atan2(reShape(:,2),reShape(:,1));
compOut = compOut.*exp(i.*phase).*(-i);
assert(angle(compOut),zeros(numDays,1),2*eps);

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
