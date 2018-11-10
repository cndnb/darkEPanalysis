pkg load signal;
%%%%%%%%%%%%%%%%%%%% PROBLEM LAYOUT & CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Damped oscillator differential equation: x'' + (2wZ) x' + (w^2) x = 0
%Where w is the undamped angular frequency, and Z = 1/2Q. Oscillation is underdamped,
%therefore solution is of the form x = e^(-wZt)[A sin(Ct) + B cos(Ct)]
%C = sqrt[w^2-(wZ)^2] = w*sqrt[1-Z^2]

%torque = kappa*theta + rotI*(d^2 theta/dt^2)

%Pendulum and balance parameters, in SI units:
global I = 378/(1e7);                                                                    
%global f0 = 1.9338e-3; %Fake data f0
global f0 = 0.0019295; %Real data f0
global Q = 500000;                                                                     
global Temp = 273+24;  
global kappa = (2*pi*f0)^2 * I;

%Sidereal day frequency
global omegaEarth = 2*pi*(1/86164.0916);

%Solar day frequency
global oED = 2*pi*(1/86400);

%Specific pendulum position data
%%Latitude, longitude, and compass direction to be input as decimal degrees
seattleLat = 47.6593743;
seattleLong = -122.30262920000001;
compassDir = 30;
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
seattleLong = (180/pi)*((pi/180)*(seattleLong + vernalEqLong)-omegaEarth*6939300);

%%%%%%%%%%%%%%%%%%%%%%%%%%% FITTER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Time that the data starts from relative to January 1st 2000 at 00:00 UTC
startTime = 0;

%Number of design matrix columns
numBETAVal = columns(createSineComponents(1,1,seattleLat,seattleLong,compassDir,startTime));
%Linear terms need constant subtracted off, need to know which column this will
linearColumn = numBETAVal - 3;

%Start of frequency scan
startFreq = 1e-4;
%End frequency scan
stopFreq = 5e-4;

%Number of days in the data considered
daysInclude = 0;

%How many periods of the specific frequency are included in weighted error fit
chunkSize = 10;

%boolean to show output on terminal
showOut = 1;

%Distance in time between samples
sampleInterval = 1; %seconds

%Boolean to use torsion filter
torsionFiltered = 1;

%Boolean for type of damping
isExternal = 0;

%Autocollimator noise floor
aCN = 1e-9;

%% Earthquake removal parameters %%

%number of seconds in a bin
dayLength = 2*86164; %seconds

%This is the value above which torque is considered an earthquake
baseThreshold = 2e-12;

%Number of seconds around a large torque that will be removed
areaRemove = [10000,20000];


%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pkg load signal;

if (!exist('d'))
	disp('importing data');
	fflush(stdout);
  importrealdarkEP
  %importfakeDarkEP
  %d = O;
endif

newD = d;
%newD = [d(1:140000,:);d(146000:2*86164,:)];

%newD = [d(1:100000,:);d(2*86164:4*86164,:)];

%newD = [d(2*86164:6*86164,:)];
%newD = d(2*86164:4*86164,:);

%newD = d(1:6*86164-20000,:); 
%newD = [d(21000:86164-20000,:);d(2*86164:3*86164-20000,:);d(3*86164:4*86164,:);d(4*86164:5*86164,:);d(5*86164+30000:6*86164-20000,:)];
%newD = d(4*86164+25000:5*86164,:);
%newD = d(1:2*86164,:);

%newD = d(1:6*86164-20000,:); 
%newD = [d(21000:86164-20000,:);d(2*86164:3*86164-20000,:);d(3*86164:4*86164,:);d(4*86164:5*86164,:);d(5*86164+30000:6*86164-20000,:)];
%newD = d(1:5*86164,:);
%newD = d;
%newD = [d(4*86164 + 22000:5*86164,:);d(5*86164+30000:6*86164-20000,:)];
%newD = blah;
%newD = [d(:,1),d(:,2)];
%newD = [d(195000:235000,:);d(370000:410000,:)];
%newD = d(370000:410000,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%% EARTHQUAKE REMOVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%returns torques set to zero at earthquakes in a matrix, 
%driftFix = data divided into days and earthquake points removed
%driftFix{day,1} = [seconds, displacement amplitude]
%Full length is length of the data in seconds from start to stop, before
%earthquake removal

noEarthquakes = 0; editTorque = 0; torqueCell = cell(2,2);
%Calculates the torque at each point, puts into an array for analysis
if(torsionFiltered)
	tFnewD = torsionFilter(newD(:,1),newD(:,2),1/f0);
	calcTorque = torque(tFnewD, I, kappa);
	threshold = baseThreshold + mean(calcTorque(3:(end-2),2));
	torqueCell{1,1} = calcTorque; torqueCell{1,2} = threshold;
	[noEarthquakes,editTorque] = removeEarthquakes(tFnewD,torqueCell,newD(:,3),areaRemove,showOut);
else
	calcTorque  = torque(newD,   I, kappa);
	threshold = baseThreshold + mean(calcTorque(3:(end-2),2));
	torqueCell{1,1} = calcTorque; torqueCell{1,2} = threshold;
	[noEarthquakes,editTorque] = removeEarthquakes(  newD,torqueCell,newD(:,3),areaRemove,showOut);
endif

%Divides data into day length chunks
driftFix  = dayDivision(noEarthquakes, daysInclude,dayLength,showOut);

removed = 0;
count = rows(driftFix);
while(count > 0)
	if(rows(driftFix{count,1}) < 86164)
		driftFix(count,:) = [];
		removed = removed + 1;
	endif
	count = count - 1;
endwhile

%Checks the total length of the data
removed
checkLength = cell2mat(driftFix(:,1));
fullLength = checkLength(end,1) - checkLength(1,1);

%Prepares data with resonance and drift removed
noRes = driftFix;
for count = 1:rows(driftFix)
	t = driftFix{count,1}(:,1);
	X = 0;
	if(torsionFiltered)
		X = [ones(rows(t),1),t-t(1,1).*ones(rows(t),1),sin(oED.*t),cos(oED.*t)];
	else
		X = [ones(rows(t),1),t-t(1,1).*ones(rows(t),1),sin(oED.*t),cos(oED.*t),sin((2*pi*f0).*t),cos((2*pi*f0).*t)];
	endif
	[bRE,sRE,rRE,errRE,covRE] = ols2(driftFix{count,1}(:,2),X);
	noRes{count,1}(:,2) = noRes{count,1}(:,2) - X*bRE;
endfor

checkRes = cell2mat(noRes);

figure(6);
plot(torqueCell{1,1}(3:(end - 2),1),torqueCell{1,1}(3:(end - 2),2),torqueCell{1,1}(3:(end - 2),1),...
torqueCell{1,2}.*ones(rows(torqueCell{1,1})-4,1));
title('Threshold plotted on torque');
xlabel('Time (s)');
ylabel('Torque (N m)');

figure(8);
plot(checkRes(:,1),checkRes(:,2));
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

[compAvg,modErr] = dispAmpTF(driftFix,freqArray,linearColumn,torsionFiltered,showOut,seattleLat,seattleLong,compassDir,startTime);

%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSION TO TORQUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sums in quadrature amplitudes to find single value for each coordinate direction,
%Divides by the transfer function to find the torque amplitude for each frequency
[FINALAMP, FINALERR,FINALPHASE] = ampToPower(compAvg,freqArray,kappa,f0,Q,sampleInterval,torsionFiltered,isExternal);

%Thermal noise limit calculation for torque and g_{B-L}
thNoise = thermalNoise(FINALAMP(:,1),kappa,Q,Temp,f0,aCN,rows(checkLength),isExternal);
gLim = torqueToGBL(thNoise);

fileName = ['run0160',ctime(time)([5:7,9:10,21:end-1]),'(',num2str(startFreq),'-',num2str(stopFreq),')'];
save(fileName,'FINALAMP','FINALERR','gLim');
disp(['Saved run to file ',fileName]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Plots torque power as a function of frequency
figure(1);
loglog(FINALAMP(:,1),[FINALAMP(:,2),thNoise]);
%%hold on;
%%loglog([FINALAMP(:,1),FINALAMP(:,1)],[FINALAMP(:,2) - modErr(:,1),FINALAMP(:,2) + modErr(:,1)],'-r');
%%hold off;
legend('Amplitude','Thermal Limit');
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency parallel to Z');

figure(2);
loglog(FINALAMP(:,1),[FINALAMP(:,3),thNoise]);
%%hold on;
%%loglog([FINALAMP(:,1),FINALAMP(:,1)],[FINALAMP(:,3) - modErr(:,2),FINALAMP(:,3) + modErr(:,2)],'-r');
%%hold off;
legend('Amplitude','Thermal Limit');
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency perpendicular to gamma');

figure(3);
loglog(FINALAMP(:,1),[FINALAMP(:,4),thNoise]);
%%hold on;
%%loglog([FINALAMP(:,1),FINALAMP(:,1)],[FINALAMP(:,4) - modErr(:,3),FINALAMP(:,4) + modErr(:,3)],'-r');
%%hold off;
legend('Amplitude','Thermal Limit');
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
