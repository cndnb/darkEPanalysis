%%%%%%%%%%%%%%%%%%%% PROBLEM LAYOUT & CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Damped oscillator differential equation: x'' + (2wZ) x' + (w^2) x = 0
%Where w is the undamped angular frequency, and Z = 1/2Q. Oscillation is underdamped,
%therefore solution is of the form x = e^(-wZt)[A sin(Ct) + B cos(Ct)]
%C = sqrt[w^2-(wZ)^2] = w*sqrt[1-Z^2]

%torque = kappa*theta + rotI*(d^2 theta/dt^2)

%Pendulum and balance parameters, in SI units:
global I = 378/(1e7);                                                                    
global f0 = 1.9338e-3; %Fake data f0
%global f0 = 0.0019295; %Real data f0
global Q = 500000;                                                                     
global Temp = 273+24;  
global kappa = (2*pi*f0)^2 * I;

%%Latitude, longitude, and compass direction to be input as decimal degrees
%global seattleLat = 47.6593743;
%global seattleLong = -122.30262920000001;
%global compassDir = 90;
%global dipoleMag = 1;
%
%%Sidereal day frequency
%global omegaEarth = 2*pi*(1/86164.0916);
%
%%Defining the X vector at January 1st 2000 00:00 UTC
%%Using the website https://www.timeanddate.com/worldclock/sunearth.html
%%At the time Mar 20 2000 07:35 UT
%%80*24*3600 + 7*3600 + 35*60 = 6939300 seconds since January 1, 2000 00:00:00 UTC
%%At the vernal equinox, longitude is equal to zero, so z=0;
%global vernalEqLat = 68.1166667;
%
%%Prepares seattleLat in terms of equatorial cordinates at January 1, 2000 00:00:00 UTC
%%This is the angle of seattleLat from the X vector
%seattleLat = rad2deg(deg2rad(seattleLat + vernalEqLat)-omegaEarth*6939300);

%%%%%%%%%%%%%%%%%%%%%%%%%%% FITTER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of design matrix columns
numBETAVal = columns(createSineComponents(1,1));
%Linear terms need constant subtracted off, need to know which column this will
%be performed on--in this analysis, it is second to last.
linearColumn = numBETAVal - 3;

%Start of frequency scan
startFreq = 1e-3;
%End frequency scan
stopFreq = 1e-2;

%Number of days in the data considered
daysInclude = 0;

%How many periods of the specific frequency are included in weighted error fit
chunkSize = 10;

%boolean to show output on terminal
showOut = 1;

%Distance in time between samples
sampleInterval = 1; %seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pkg load signal;

if (!exist('d'))
  %importrealdarkEP
  %importfakeDarkEP
  d = O;
endif


%newD = [d(1:140000,:);d(146000:2*86164,:)];

newD = d(1:2*86164,:);

%newD = d(1:6*86164-20000,:); 
%newD = [d(21000:86164-20000,:);d(2*86164:3*86164-20000,:);d(3*86164:4*86164,:);d(4*86164:5*86164,:);d(5*86164+30000:6*86164-20000,:)];
%newD = d(1:5*86164,:);
%newD = d;
%newD = [d(4*86164 + 22000:5*86164,:);d(5*86164+30000:6*86164-20000,:)];
%newD = blah;
%newD = [d(:,1),d(:,2)];
%newD = [d(195000:235000,:);d(370000:410000,:)];
%newD = d(370000:410000,:);

tFnewD = torsionFilter(newD(:,1),newD(:,2),1/f0);
%%%%%%%%%%%%%%%%%%%%%%%%%%% EARTHQUAKE REMOVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculates the torque at each point, puts into an array for analysis
calcTorque  = torque(newD,   I, kappa);
calcTorque2 = torque(tFnewD, I, kappa);


%number of seconds in a bin
dayLength = 86164; %seconds
%This is the value above which torque is considered an earthquake
threshold = 1e-12 + mean(calcTorque(3:(end-2),2));
%threshold = 1e-13 + mean(calcTorque(3:(end-2),2));
%Number of seconds around a large torque that will be removed

areaRemove = [11000,25000];


%returns torques set to zero at earthquakes in a matrix, 
%driftFix = data divided into days and earthquake points removed
%driftFix{day,1} = [seconds, displacement amplitude]
%Full length is length of the data in seconds from start to stop, before
%earthquake removal

[noEarthquakes,editTorque]   = removeEarthquakes(newD,  calcTorque, threshold,areaRemove,showOut);
[noEarthquakes2,editTorque2] = removeEarthquakes(tFnewD,calcTorque2,threshold,areaRemove,showOut);


omegaEarth = 2*pi*(1/86164.0916);
oED = 2*pi*(1/86400);
t = noEarthquakes(:,1);
X = [ones(rows(t),1),t-t(1,1).*ones(rows(t),1),sin(oED.*t),cos(oED.*t)];
[bRE,sRE,rRE,errRE,covRE] = ols2(noEarthquakes(:,2),X);
noRes = [noEarthquakes(:,1),noEarthquakes(:,2) - X*bRE];

driftFix  = dayDivision(noEarthquakes, daysInclude,dayLength,showOut);
driftFix2 = dayDivision(noEarthquakes2,daysInclude,dayLength,showOut);

checkLength = cell2mat(driftFix(:,1));
fullLength = checkLength(end,1) - checkLength(1,1);

checkLength2 = cell2mat(driftFix2(:,1));
fullLength2 = checkLength2(end,1) - checkLength2(1,1);


figure(6);
plot(calcTorque(3:(end - 2),1),calcTorque(3:(end - 2),2),calcTorque(3:(end - 2),1),...
threshold.*ones(rows(calcTorque)-4,1),editTorque(:,1),editTorque(:,2));
title('Threshold plotted on torque');
xlabel('Time (s)');
ylabel('Torque (N m)');

figure(7);
plot(noRes(:,1),noRes(:,2));
title('drift, earthRotation, and resonance removed');
xlabel('Time (s)');
ylabel('Angle');

%%%%%%%%%%%% AMPLITUDE(TIME) => AMPLITUDE(FREQUENCY) CONVERSION  %%%%%%%%%%%%%%%%

%Fits design matrix at each frequency over a given number of data bins
%ampFreq is the average of each amplitude at each frequency over all bins
%ampError is the standard deviation of each amplitude at each frequency over all bins

%endCount is the #rows of frequency matrix
%endCout = total frequency band divided by the smallest frequency jump
%Integer so that it can be used for indexing

freqArray = 1:(floor(fullLength/2));
freqArray = [0,freqArray];
freqArray = freqArray';
freqArray = freqArray./fullLength;

tempStart = freqArray - startFreq.*ones(rows(freqArray),1);
tempEnd = freqArray - stopFreq.*ones(rows(freqArray),1);
pastMinStart = Inf;
pastMinEnd = Inf;
minIndStart = 0;
minIndEnd = 0;
for count = 1:rows(freqArray)
  if (abs(tempStart(count)) < pastMinStart)
    pastMinStart = abs(tempStart(count));
    minIndStart = count;
  endif
  if (abs(tempEnd(count)) < pastMinEnd)
    pastMinEnd = abs(tempEnd(count));
    minIndEnd = count;
  endif
endfor
indStart = minIndStart;
indEnd = minIndEnd;
freqArray = freqArray(indStart:indEnd,1);

freqArray(1,1)
freqArray(end,1)
rows(freqArray)
fflush(stdout);
pause();

[compAvg,compOut] = dispAmpTF(driftFix,freqArray,linearColumn,0,showOut);


%Second interation
freqArray2 = 1:(floor(fullLength2/2));
freqArray2 = [0,freqArray2];
freqArray2 = freqArray2';
freqArray2 = freqArray2./fullLength2;

tempStart = freqArray2 - startFreq.*ones(rows(freqArray2),1);
tempEnd = freqArray2 - stopFreq.*ones(rows(freqArray2),1);
pastMinStart = Inf;
pastMinEnd = Inf;
minIndStart = 0;
minIndEnd = 0;
for count = 1:rows(freqArray2)
  if (abs(tempStart(count)) < pastMinStart)
    pastMinStart = abs(tempStart(count));
    minIndStart = count;
  endif
  if (abs(tempEnd(count)) < pastMinEnd)
    pastMinEnd = abs(tempEnd(count));
    minIndEnd = count;
  endif
endfor
indStart = minIndStart;
indEnd = minIndEnd;
freqArray2 = freqArray2(indStart:indEnd,1);

freqArray2(1,1)
freqArray2(end,1)
rows(freqArray2)
fflush(stdout);
pause();

[compAvg2,compOut2] = dispAmpTF(driftFix2,freqArray2,linearColumn,1,showOut);

%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSION TO TORQUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sums in quadrature amplitudes to find single value for each coordinate direction,
%Divides by the transfer function to find the torque amplitude for each frequency
[FINALAMP, FINALERR,FINALPHASE] = ampToPower(compAvg,freqArray,kappa,f0,Q);
[FINALAMP2, FINALERR2,FINALPHASE2] = ampToPower(compAvg2./twoPointTransfer(freqArray2,f0,1),freqArray2,kappa,f0,Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plots torque power as a function of frequency
figure(1);
loglog(FINALAMP(:,1),FINALAMP(:,2),FINALAMP2(:,1),FINALAMP2(:,2));
%%hold on;
%%loglog([FINALAMP(:,1),FINALAMP(:,1)],[FINALAMP(:,2) - FINALERR(:,2),FINALAMP(:,2) + FINALERR(:,2)],'-r');
%%hold off;
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency parallel to Z');

figure(2);
loglog(FINALAMP(:,1),FINALAMP(:,3),FINALAMP2(:,1),FINALAMP2(:,3));
%%hold on;
%%loglog([FINALAMP(:,1),FINALAMP(:,1)],[FINALAMP(:,3) - FINALERR(:,3),FINALAMP(:,3) + FINALERR(:,3)],'-r');
%%hold off;
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency perpendicular to gamma');

figure(3);
loglog(FINALAMP(:,1),FINALAMP(:,4),FINALAMP2(:,1),FINALAMP2(:,4));
%%hold on;
%%loglog([FINALAMP(:,1),FINALAMP(:,1)],[FINALAMP(:,4) - FINALERR(:,4),FINALAMP(:,4) + FINALERR(:,4)],'-r');
%%hold off;
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency parallel to gamma');

%figure(4);
%loglog(FINALAMP(:,1),[FINALAMP(:,2),FINALAMP(:,3),FINALAMP(:,4)]);
%legend('Z component','Perpendicular to gamma','Parallel to gamma');
%xlabel('Frequency (Hz)');
%ylabel('Torque (N m)');
%title('Torque vs frequency');


%!test
%! disp('removeEarthquakes');
%! test removeEarthquakes
%! disp('dayDivision')
%! test dayDivision
%! disp('dispAmpTF');
%! test dispAmpTF
%! disp('createSineComponents');
%! test createSineComponents
%! disp('frequencyVariance');
%! test frequencyVariance
%! disp('specFreqAmp');
%! test specFreqAmp
%! disp('weightedOLS');
%! test weightedOLS
%! disp('transferFunction');
%! test transferFunction
%! disp('twoPointTransfer');
%! test twoPointTransfer
%! disp('ampToPower');
%! test ampToPower
%! disp('testerFile');
%! fflush(stdout);
