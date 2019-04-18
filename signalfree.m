%'torqueSim'
%test torqueSim


startTime = 0;
global omegaEarth = 2*pi*(1/86164.0916);
%Specific pendulum position data
%%Latitude, longitude, and compass direction to be input as decimal degrees
degSeattleLat = 47.6593743;
degSeattleLong = -122.30262920000001;
degCompassDir = 30;
%global dipoleMag = 1;
%
%%Defining the X vector at January 1st 2000 00:00 UTC
%%Using the website https://www.timeanddate.com/worldclock/sunearth.html
%%At the time Mar 20 2000 07:35 UT
%%80*24*3600 + 7*3600 + 35*60 = 6939300 seconds since January 1, 2000 00:00:00 UTC
%%At the vernal equinox, longitude is equal to zero, so z=0;
global vernalEqLong = 68.1166667;
%
%%Prepares seattleLat in terms of equatorial cordinates at January 1, 2000 00:00:00 UTC
%%This is the angle of seattleLat from the X vector
%%Additionally converts all angles to radians
seattleLong = ((pi/180)*(degSeattleLong + vernalEqLong)-omegaEarth*6939300);
seattleLat = (pi/180)*degSeattleLat;
compassDir = (pi/180)*degCompassDir;

%Make some time
t = 1:2*86400;t=t';


A = 1e-18;
f = 5e-3;
%omegaEarth = 2*pi*(1/86164.0916);

cM = preCalcComponents(t,seattleLat,seattleLong,compassDir,startTime);
cS = [1,1,0,0,0,0,0,0,0,0];
cS2 = [0,0,1,1,0,0,0,0,0,0];

finalSignal = sum(A*createSineComponents(t,f,cM,cS),2)+ sum(A*createSineComponents(t,7/5*f,cM,cS2));

%Parameters of the experiment
I = 378/1e7;                                                                    
f0 = 1.9338e-3;                                                                 
kappa = ((2*pi*f0)^2)*I;                                                        
Q = 500000;                                                                     
Temp = 0; 

%Simulate a pendulum
T= torqueSim(t,I, kappa,Q, Temp, finalSignal);

%Fake some autocollimator noise
AutocollimatorNoise = randn(size(t)) * 0.5e-9;

%Generate measured angle output (THIS IS THE SAME AS THE FAKE DATASET!)
O = [T(:,1) T(:,2) ];%+ AutocollimatorNoise];

%Re-compute torque
accel = diff(diff(O(:,2)));
Tor = I*accel + kappa*O(2:end-1,2);

fullLength = rows(O);

%Checks that peaks are at the correct points
figure(1);
check = psd(t(2:length(t)-1,1),Tor);
loglog(check(:,1),check(:,2));


