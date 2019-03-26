%'torqueSim'
%test torqueSim

%Make some time
t = 1:2*86400;t=t';


A = 1e-16;
f = 9e-3;
%omegaEarth = 2*pi*(1/86164.0916);

cM = preCalcComponents(t,seattleLat,seattleLong,compassDir,startTime);
cS = [1,0,0,0,0,0,0,0,0,0];

finalSignal = A*createSineComponents(t,f,cM,cS);
assert(columns(finalSignal) == 1);

%Parameters of the experiment
I = 378/1e7;                                                                    
f0 = 1.9338e-3;                                                                 
kappa = ((2*pi*f0)^2)*I;                                                        
Q = 500000;                                                                     
Temp = 273+24; 

%Simulate a pendulum
T= torqueSim(t,I, kappa,Q, Temp, finalSignal);

%Fake some autocollimator noise
AutocollimatorNoise = randn(size(t)) * 0.5e-9;

%Generate measured angle output (THIS IS THE SAME AS THE FAKE DATASET!)
O = [T(:,1) T(:,2) + AutocollimatorNoise];

%Re-compute torque
accel = diff(diff(O(:,2)));
Tor = I*accel + kappa*O(2:end-1,2);

fullLength = rows(O);

%Checks that peaks are at the correct points
figure(1);
check = psd(t(2:length(t)-1,1),Tor);
loglog(check(:,1),check(:,2));


