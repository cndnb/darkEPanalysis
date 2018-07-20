
%t=1:1e5; t=t';
%tAmp = 1e-15;
%signal = tAmp.*sin((2*pi*1e-3).*t);
%
%data = [t,signal];

data = O(1:100000,:);
I = 378/(1e7);                                                                    
f0 = 1.9338e-3;                                                                 
Q = 500000;                                                                     
T = 273+24;  
kappa = ((2*pi*f0)^2) * I;



ampFreq = ones(length(5e-4:(1/rows(data)):5e-3),2);
ampFreqOLS = ampFreq;
count = 1;

BETAARR = ones(rows(ampFreq),1);
BETAOLSARR = ones(rows(ampFreq),1);
for freq = 1e-3:(1/rows(data)):1e-2
  freq
  fflush(stdout);
[BETA,COV] = specFreqAmp(data,createSineComponents(data(:,1),freq),sW);
[BETAOLS,COVOLS] = ols2(data(:,2),createSineComponents(data(:,1),freq));
BETAOLS = BETAOLS';
ampFreq(count,:) = [freq,sqrt(BETA(1)^2+BETA(2)^2)];
ampFreqOLS(count,:) = [freq,sqrt(BETAOLS(1)^2+BETAOLS(2)^2)];
ampFreq(count,2) = abs(ampFreq(count,2)/transferFunction(freq,kappa,f0,Q));
ampFreqOLS(count,2) = abs(ampFreqOLS(count,2)/transferFunction(freq,kappa,f0,Q));
count = count+1;
endfor

figure(1);
plot(data(:,1),data(:,2));

figure(2);
loglog(ampFreq(:,1),ampFreq(:,2),ampFreqOLS(:,1),ampFreqOLS(:,2));
legend('specFreqPower','OLS2');

figure(3);
check = psd(data(:,1),data(:,2));
loglog(check(:,1),check(:,2));