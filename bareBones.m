%if(!exist('d'))
%	importfakeDarkEP
%endif
%cData = d(370001:410000,:);
if(!exist('O'))
  signalfree
endif
cData = O;

seattleLat = pi/4;
seattleLong = pi/4;
compassDir = pi/6;
startTime = 0;

I = 378/(1e7);                                                                    
f0 = 1.9338e-3;                                                                 
Q = 500000;                                                                     
Temp = 273+24;  
kappa = (2*pi*f0)^2 * I;
omegaEarth = 2*pi*(1/86164.0916);

freqArray = (0:(rows(cData)/2))/(rows(cData));
freqArray = freqArray';
tau = transferFunction(freqArray,kappa,f0,Q,0);

dFX = [ones(rows(cData),1),(1:rows(cData))',sin(cData(:,1).*omegaEarth),cos(cData(:,1).*omegaEarth)];
[dFB,dFS,dFR,dFERR,dFCOV] = ols2(cData(:,2),dFX);
cData(:,2) = cData(:,2) - dFX*dFB;

%Finds FFT output
freqData = (2/rows(cData))*fft(cData(:,2));
assert(freqData(2:rows(freqData)/2+1,:),conj(flip(freqData((rows(freqData)/2)+1:end,:))));
freqData = freqData(1:(rows(freqData)/2+1),:);
torqueData = abs(freqData./tau);


startVal = 99;
freqArray(startVal+2)
endVal = 48999;
freqArray(rows(freqArray) - 1 - endVal)
rows(freqArray) - endVal - 1
fflush(stdout);
pause();

collectionArray = ones(rows(freqArray),3);
cAS = ones(rows(freqArray),3);
covArray = ones(12,12,rows(freqArray));
rows(freqArray)
for freq = 2+startVal:rows(freqArray) - 1 - endVal
	freq
	fflush(stdout);
	designX = createSineComponents(cData(:,1),freqArray(freq),seattleLat,seattleLong,compassDir,startTime);
	[b,s,r,err,cov] = ols2(cData(:,2),designX);
	covArray(:,:,freq) = cov;
	[bZ,s,r,err,cov] = ols2(cData(:,2),designX(:,1:2));
	[bPeX,s,r,err,cov] = ols2(cData(:,2),designX(:,3:4));
	[bPaX,s,r,err,cov] = ols2(cData(:,2),designX(:,5:6));
	collectionArray(freq,:) = [b(2,1) + i.*b(1,1),b(4,1) + i.*b(3,1),b(6,1) + i.*b(5,1)];
	cAS(freq,:) = [bZ(2,1) + i.*bZ(1,1),bPeX(2,1) + i.*bPeX(1,1),bPaX(2,1) + i.*bPaX(1,1)];
endfor
ampOut = [freqArray,abs(collectionArray./tau)];
aOS = [freqArray,abs(cAS./tau)];
figure(1);
title('Joint fit');
loglog(ampOut(:,1),[ampOut(:,2),ampOut(:,3),ampOut(:,4),torqueData]);
legend('Z','PerpX','ParaX','FFT');
figure(2);
title('Separate fits');
loglog(aOS(:,1),[aOS(:,2),aOS(:,3),aOS(:,4),torqueData]);
legend('Z','PerpX','ParaX','FFT');

