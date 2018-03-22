d = load('fakeDarkEPAugust92017.dat');

I = 378/(1e7);                                                                    
f0 = 1.9338e-3;                                                                 
Q = 500000;                                                                     
T = 273+24;  
kappa = (2*pi*f0)^2 * I;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRIFT CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Removes drift from the data
linearDesign = [ones(length(d(:,1)),1),d(:,1)];
[lDB, lDS, lDR, LDEr, lDCov] = ols2(d(:,2),linearDesign);
oldDriftFix = [d(:,1),d(:,2)-(linearDesign*lDB)];

%Removes the daily swing from the data
[dailyBeta, dailySigma, dailyR] = sineFitter(oldDriftFix(:,1),oldDriftFix(:,2),1/86400);
oldDriftFix = [oldDriftFix(:,1),oldDriftFix(:,2)-genSineSeed(oldDriftFix(:,1),1/86400)*dailyBeta];

%%%%%%%%%%%%%%%%%%%%%%%% EARTHQUAKE REMOVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This is the value above which torque is considered an earthquake
threshold = 1e-13;
areaRemove = 10000;

calcTorque = torque(oldDriftFix, I, kappa);
[driftFix, editTorque,pointsRemoved] = removeEarthquakes(oldDriftFix,calcTorque,threshold,areaRemove);
check = psd(editTorque(:,1),editTorque(:,2));



figure(1);
plot(oldDriftFix(:,1),oldDriftFix(:,2),driftFix(:,1),driftFix(:,2)+2e-5);

%%Checking that FFT of torque has no peaks
%figure(1);
%loglog(check(:,1),check(:,2));
%
%Checking the threshold level
figure(2);
plot(calcTorque(3:(1e6 - 2),1),calcTorque(3:(1e6 - 2),2),calcTorque(3:(1e6 - 2),1),...
threshold.*ones(length(calcTorque)-4,1));