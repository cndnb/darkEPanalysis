dataDivisions = cell(2,1);
f = 2*pi*(5e-3);
dataDivisions{1,1} = [(1:40000)',sin(f.*(1:40000)')+randn(40000,1).*10];
dataDivisions{2,1} = [(105001:145000)',sin(f.*(105001:145000)')+randn(40000,1).*20];
startFreq = 1e-3;
stopFreq = 1e-2;

fullLength = 40000;

endCount = floor((stopFreq-startFreq)/(jump*(1/fullLength)))+1
fflush(stdout);

freqArray = ones(endCount,1);
  
%Assigns frequency values for the first column of the frequency and error arrays
for count = 1:endCount
  freqArray(count,1) = (startFreq+((count-1)*(1/fullLength))); %fullLength passed before earthquakes removed
endfor

nDF = driftFix;
nDF{1,1} = nDF{1,1}(10000:end,:);

[ampFreq1,ampError1] = dispAmpTF(nDF(1,1),freqArray,endCount,0,0,1);
[ampFreq2,ampError2] = dispAmpTF(nDF(2,1),freqArray,endCount,0,0,1);

pwr = ones(rows(ampFreq1),2);
pwr(:,1) = sqrt(ampFreq1(:,1).^2 + ampFreq1(:,2).^2);
pwr(:,2) = sqrt(ampFreq2(:,1).^2 + ampFreq2(:,2).^2);
pwr(:,3) = (pwr(:,1) + pwr(:,2))./2;

figure(1);
loglog(freqArray,pwr(:,1),freqArray,pwr(:,2),freqArray,pwr(:,3))
figure(2);
loglog(freqArray,abs(pwr(:,3)./transferFunction(freqArray,kappa,f0,Q)))