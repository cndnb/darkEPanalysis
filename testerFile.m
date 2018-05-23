%%%%%%%%%%%%%%%%%%%% UNIT TESTS/PRELIMINARY CHECKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'dispAmpTF'
test dispAmpTF
'createSineComponents'
test createSineComponents
'frequencyVariance'
test frequencyVariance
'specFreqPower'
test specFreqAmp
'weightedOLS'
test weightedOLS
'transferFunction'
test transferFunction
'ampToPower'
test ampToPower

fflush(stdout);
pause()

%%%%%%%%%%%%%%%%%%%% PROBLEM LAYOUT & CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Damped oscillator differential equation: x'' + (2wZ) x' + (w^2) x = 0
%Where w is the undamped angular frequency, and Z = 1/2Q. Oscillation is underdamped,
%therefore solution is of the form x = e^(-wZt)[A sin(Ct) + B cos(Ct)]
%C = sqrt[w^2-(wZ)^2] = w*sqrt[1-Z^2]

%torque = kappa*theta + rotI*(d^2 theta/dt^2)

%Pendulum and balance parameters, in SI units:
I = 378/(1e7);                                                                    
f0 = 1.9338e-3;                                                                 
Q = 500000;                                                                     
T = 273+24;  
kappa = (2*pi*f0)^2 * I;

%Latitude, longitude, and compass direction to be input as decimal degrees
global seattleLat = 47.6593743;
global seattleLong = -122.30262920000001;
global compassDir = 90;
global dipoleMag = 1;

%Sidereal day frequency
global omegaEarth = 2*pi*(1/86164.0916);

%Defining the X vector at January 1st 2000 00:00 UTC
%Using the website https://www.timeanddate.com/worldclock/sunearth.html
%At the time Mar 20 2000 07:35 UT
%80*24*3600 + 7*3600 + 35*60 = 6939300 seconds since January 1, 2000 00:00:00 UTC
%At the vernal equinox, longitude is equal to zero, so z=0;
global vernalEqLat = 68.1166667;

%Prepares seattleLat in terms of equatorial cordinates at January 1, 2000 00:00:00 UTC
%This is the angle of seattleLat from the X vector
seattleLat = rad2deg(deg2rad(seattleLat + vernalEqLat)-omegaEarth*6939300);

%%%%%%%%%%%%%%%%%%%%%%%%%%% FITTER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1 for weightedOLS, 0 for ols2
fitIsWeighted = 1;

%Number of design matrix columns
numBETAVal = 12;
%How many periods of the specific frequency are included in error fit
chunkSize = 50;
%Multiples of smallest usable frequency between amplitude points
jump = 1000;
%Start of frequency scan
startFreq = 1e-3;
%End frequency scan
stopFreq = 1e-2;
%How many points are included in the coherent average bins
dataDivisions = 5;

%%%%%%%%%%%% TIME => FREQUENCY DEPENDENT DISPLACEMENT AMPLITUDE %%%%%%%%%%%%%%%%

%Fits design matrix at each frequency over a given number of data bins
%ampFreq is the average of each amplitude at each frequency over all bins
%ampError is the standard deviation of each amplitude at each frequency over all bins

%endCount is the #rows of frequency matrix
%endCout = total frequency band divided by the smallest frequency jump
%Integer so that it can be used for indexing
endCount = floor((stopFreq-startFreq)/(jump*(1/fullLength)))+1;

freqArray = ones(endCount,1);
  
%Assigns frequency values for the first column of the frequency and error arrays
for count = 1:endCount
  freqArray(count,1) = (startFreq+((count-1)*jump*(1/fullLength))); %fullLength passed before earthquakes removed
endfor
  
[ampFreq,ampError] = dispAmpTF(driftFix,freqArray,endCount,dataDivisions,chunkSize,numBETAVal,0,fitIsWeighted);

%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSION TO TORQUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creates array for final data
FINALAMP = ones(rows(ampFreq),5);
FINALERR = ones(rows(ampError),4);

%Sums in quadrature amplitudes to find single value for each coordinate direction,
%Divides by the transfer function to find the torque amplitude for each frequency
[FINALAMP, FINALERR] = ampToPower(ampFreq,ampError,kappa,f0,Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plots torque power as a function of frequency
figure(5);
loglogerr(FINALAMP(:,1),FINALAMP(:,2),FINALERR(:,2));
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency parallel to gamma');

figure(6);
loglogerr(FINALAMP(:,1),FINALAMP(:,3),FINALERR(:,3));
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency perpendicular to gamma');

figure(7);
loglogerr(FINALAMP(:,1),FINALAMP(:,4),FINALERR(:,4));
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency in the z component');

figure(8);
loglog(FINALAMP(:,1),FINALAMP(:,2),FINALAMP(:,1),FINALAMP(:,3),FINALAMP(:,1),FINALAMP(:,4),FINALAMP(:,1),FINALAMP(:,5));
legend('Parallel to gamma','Perpendicular to gamma','Z component','Sum signal');
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency');