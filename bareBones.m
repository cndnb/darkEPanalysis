importfakeDarkEP
cData = d(370000:410000,:);
freqArray = (1:(rows(cData)-1))/(2*rows(cData));
freqArray = freqArray';

dFX = [ones(rows(cData),1),cData(:,1),sin(cData(:,1).*omegaEarth),cos(cData(:,1).*omegaEarth)];
[dFB,dFS,dFR,dFERR,dFCOV] = ols2(cData(:,2),dFX);
cData(:,2) = cData(:,2) - dFX*dFB;

collectionArray = ones(rows(freqArray),3);
rows(freqArray)
for freq = 1:rows(freqArray)
	freq
	fflush(stdout);
	designX = createSineComponents(cData(:,1),freqArray(freq));
	[b,s,r,err,cov] = ols2(cData(:,2),designX);
	collectionArray(count,:) = [b(2,1) + i.*b(1,1),b(4,1) + i.*b(3,1),b(6,1) + i.*b(5,1)];
endfor
ampOut = [freqArray,abs(collectionArray./transferFunction(freqArray,kappa,f0,Q))];
figure(1);
loglog(ampOut(:,1),ampOut(:,2),ampOut(:,3),ampOut(:,4));

