
 b=1:1000;
 b=b';
 b=[b,randn(1000,1)];
 fitFreq = 1e-2;
 [B,olsSigma,R,Err,Cov] = ols2(b(:,2),createSineComponents(b(:,1),fitFreq));
 varRes = frequencyVariance(b,fitFreq,rows(b)*fitFreq);
 fVSigma = varRes(1);
 assert (fVSigma == (1./olsSigma))