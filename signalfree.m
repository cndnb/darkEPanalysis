%Make some time
t = 1:1e6;t=t';
searchFreq = [1e-2,1e-3,1e-4,1e-5];
f = 1e-3;
A = 1e-18;

finalSignal = A*sin(2*pi*searchFreq(1)*t)+A*sin(2*pi*searchFreq(2)*t)+...
A*sin(2*pi*searchFreq(3)*t)+A*sin(2*pi*searchFreq(4)*t);


%Parameters of the experiment
I = 378/1e7;                                                                    
f0 = 1.9338e-3;                                                                 
kappa = (2*pi)^2*f0^2*I;                                                        
Q = 500000;                                                                     
T = 273+24; 

%Simulate a pendulum
T= torqueSim(t,I, kappa,Q, T, finalSignal);

%Fake some autocollimator noise
AutocollimatorNoise = randn(size(t)) * 0.5e-9;

%Generate measured angle output (THIS IS THE SAME AS THE FAKE DATASET!)
O = [T(:,1) T(:,2) + AutocollimatorNoise];

%Re-compute torque
accel = diff(diff(O(:,2)));
Tor = I*accel + kappa*O(2:end-1,2);

%Checks that peaks are at the correct points
check = psd(t(2:length(t)-1,1),Tor);
loglog(check(:,1),check(:,2));

