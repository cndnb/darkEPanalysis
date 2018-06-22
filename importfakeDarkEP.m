d = load('fakeDarkEPAugust92017.dat');

I = 378/(1e7);                                                                    
f0 = 1.9338e-3;                                                                 
Q = 500000;                                                                     
T = 273+24;  
kappa = (2*pi*f0)^2 * I;


%%%%%%%%%%%%%%%%%%%%%%%% DRIFT CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Removes drift from the data
%linearDesign = [ones(length(d(:,1)),1),d(:,1)];
%[lDB, lDS, lDR, LDEr, lDCov] = ols2(d(:,2),linearDesign);
%oldDriftFix = [d(:,1),d(:,2)-(linearDesign*lDB)];
%
%%Removes the daily swing from the data
%[dailyBeta, dailySigma, dailyR] = sineFitter(oldDriftFix(:,1),oldDriftFix(:,2),1/86400);
%oldDriftFix = [oldDriftFix(:,1),oldDriftFix(:,2)-genSineSeed(oldDriftFix(:,1),1/86400)*dailyBeta];


%%%%%%%%%%%%%%%%%%%%%%%% EARTHQUAKE REMOVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calcTorque = torque(d, I, kappa);

%This is the value above which torque is considered an earthquake
threshold = 1e-13 + mean(calcTorque(3:(1e6-2),2));
areaRemove = 10000;

[driftFix,editTorque] = removeEarthquakes(d,calcTorque,threshold,areaRemove);
check = psd(editTorque(:,1),editTorque(:,2)-mean(editTorque(:,2)));


numRows = 0;
for count = 1:rows(driftFix)
  numRows = rows(driftFix{count,1});
endfor

fullData = zeros(numRows,2);
indexNum = 1
for count = 1:rows(driftFix)
  fullData(indexNum:(indexNum+rows(driftFix{count,1})-1),:) = driftFix{count,1};
  indexNum = indexNum + rows(driftFix{count,1});
endfor

fullLength = rows(fullData);

figure(1);
plot(d(:,1),d(:,2));
figure(2);
plot(fullData(:,1),fullData(:,2));

%Checking that FFT of torque has no peaks
figure(5);
loglog(check(:,1),check(:,2));

%Checking the threshold level
figure(6);
plot(calcTorque(3:(1e6 - 2),1),calcTorque(3:(1e6 - 2),2),calcTorque(3:(1e6 - 2),1),...
threshold.*ones(length(calcTorque)-4,1));