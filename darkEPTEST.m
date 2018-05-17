function darkEPTEST(tLength,nChunks)
t= 1:tLength; t=t';


tauT = [t,randn(length(t),1)]; %Time domain torque (random noise)

nChunks = 100;

%Pendulum constants
I = 378/1e7;                                                                    
f0 = 1.9338e-3;                                                                 
kappa = ((2*pi*f0)^2)*I;                                                        
Q = 1;                                                                     
T = 0; %No temp noise

chunkLength = floor(length(t)/nChunks);
    
preNumFreq = psd(tauT(1:1+chunkLength,:));
numFreq = rows(preNumFreq(:,1));
allTauW = ones(numFreq,1+nChunks);
allThetaW = ones(numFreq,1+nChunks);


for count = 0:nChunks-2

  preTheta = torqueSim(t(1+count*chunkLength:1+(count+1)*chunkLength,1),I,kappa,Q,T,tauT(1+count*chunkLength:1+(count+1)*chunkLength,2));

  thetaT = [preTheta(:,1),preTheta(:,2)]; %Pulls out displacement values from torqueSim

  tauW = psd(tauT(1+count*chunkLength:1+(count+1)*chunkLength,1),tauT(1+count*chunkLength:1+(count+1)*chunkLength,2)); %Frequency power spectrum of tau
  thetaW = psd(thetaT(:,1),thetaT(:,2)); %Frequency power specturm of theta
  allTauW(:,1) = tauW(:,1);
  allThetaW(:,1) = thetaW(:,1);
  allTauW(:,1+count) = tauW(:,2);
  allThetaW(:,1+count) = thetaW(:,2);
endfor

avgTauW = [allTauW(:,1),mean(allTauW(:,2:1+nChunks)')'];
avgThetaW = [allThetaW(:,1),mean(allThetaW(:,2:1+nChunks)')'];

%figure(1);
%loglog(avgTauW(:,1),avgTauW(:,2),avgThetaW(:,1),avgThetaW(:,2)); %Frequency domain comparison
% 
%figure(2);
%tfVal = abs(transferFunction(avgTauW(:,1),kappa,f0,Q).^2); %ratio of power spectra vs power of transfer function
%loglog(avgTauW(:,1),avgThetaW(:,2)./avgTauW(:,2),avgTauW(:,1),tfVal);
%
%figure(3);
%semilogx(avgTauW(:,1),avgThetaW(:,2)./(avgTauW(:,2).*tfVal)-1);

for num = 1:rows(avgTauW)-1000
  assert (abs(avgThetaW(num,2)/(avgTauW(num,2).*tfVal(num))-1)<0.07) %0.07 is the value that this data doesn't pass
endfor