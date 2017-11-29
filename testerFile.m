chunkSize = 50;
jump = 100;
startFreq = 1e-3;
stopFreq = 50e-3;
dataCut = 100000;
ampFreq = zeros((stopFreq-startFreq)/(jump*(1/dataCut)),7);
ampError = zeros((stopFreq-startFreq)/(jump*(1/dataCut)),7);
count = 0;
for j=1:(dataCut):length(driftFix)-dataCut
%Sums each data cut into one array
[sAmp,sVar] = fakeDarkEPanalysis(driftFix(j:(j+dataCut),:),chunkSize,jump,startFreq,stopFreq);
ampFreq = ampFreq + sAmp;
ampError = ampError + 1./sVar;
count = count + 1;
endfor
%Divides by the number of sums for the average
ampFreq(:,2:7) = ampFreq(:,2:7)./count;
ampError = sqrt(1./ampError);

%Creates array for final data
FINALAMP = ones(length(ampFreq),2);
FINALERR = ones(length(ampError),2);
for count = 1:length(ampFreq)
%Sums in quadrature the averaged sine/cosine amplitudes (inner product)
FINALAMP(count,:) = [ampFreq(count,1),sqrt(ampFreq(count,2:7)*ampFreq(count,2:7)')];
endfor
%Plots power over frequency
figure(1);
loglog(FINALAMP(:,1),FINALAMP(:,2));