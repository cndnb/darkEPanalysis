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

%%%%%%%%%%%%%%%%%%%%%%%%%%% FITTER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Time that the data starts from relative to January 1st 2000 at 00:00 UTC
startTime = 0;

%Number of design matrix columns
%columnSelector = [sinZ,cosZ,sinperpX,cosparaX,sinparaX,cosparaX,drift,constant,sinf0,cosf0]
%columnSelector entries = 1 will be included in the fit.
columnSelector = [1,1,1,1,1,1,1,1,0,0];

%Boolean value, 1 will run a torsion filter at the resonance frequency.
torsionFiltered = 1;

%Start of frequency scan
startFreq = 1e-3;
%End frequency scan
stopFreq = 1e-2;

%Boolean to save analyzed data to file after run is complete
saveData = 0;

%Number of days in the data considered
daysInclude = 0;

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
dayLength = 86164; %seconds

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

%This variable can be changed to remove bad data
newD = d;

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

%Gets rid of chunks with less than half of dayLength of data
removed = 0;
count = rows(driftFix);
while(count > 0)
	if(rows(driftFix{count,1}) < dayLength/2)
		driftFix(count,:) = [];
		removed = removed + 1;
	endif
	count = count - 1;
endwhile

%Checks the total length of the data 
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

%figure(6);
%plot(torqueCell{1,1}(3:(end - 2),1),torqueCell{1,1}(3:(end - 2),2),torqueCell{1,1}(3:(end - 2),1),...
%torqueCell{1,2}.*ones(rows(torqueCell{1,1})-4,1));
%title('Threshold plotted on torque');
%xlabel('Time (s)');
%ylabel('Torque (N m)');

%figure(8);
%plot(checkRes(:,1),checkRes(:,2));
%title('drift, earthRotation, and resonance removed');
%xlabel('Time (s)');
%ylabel('Angle');

%%%%%%%%%%%% AMPLITUDE(TIME) => AMPLITUDE(FREQUENCY) CONVERSION  %%%%%%%%%%%%%%%%

%Fits design matrix at each frequency over a given number of data bins
%ampFreq is the average of each amplitude at each frequency over all bins
%ampError is the standard deviation of each amplitude at each frequency over all bins

%endCount is the #rows of frequency matrix
%endCout = total frequency band divided by the smallest frequency jump
%Integer so that it can be used for indexing

%Finds index of start frequency and end frequency and modifies freqArray to run between them
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

disp("Starting Frequency:");
freqArray(1,1)
disp("Ending Frequency:");
freqArray(end,1)
disp("Number of frequencies to scan:");
rows(freqArray)
disp("Number of days");
rows(driftFix)
disp("Number of days removed");
removed
fflush(stdout);

%Pause added to frequencies can be checked, and the length of run time can be estimated. Also the torque plots should be considered to make sure that the correct data is being analyzed.
pause(5);


[compAvg,modErr] = dispAmpTF(driftFix,freqArray,columnSelector,showOut,seattleLat,seattleLong,compassDir,startTime);

%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSION TO TORQUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sums in quadrature amplitudes to find single value for each coordinate direction,
%Divides by the transfer function to find the torque amplitude for each frequency
[FINALAMP, FINALERR,FINALPHASE] = ampToPower(compAvg,modErr,freqArray,kappa,f0,Q,sampleInterval,torsionFiltered,isExternal);


%Initializes array for 95% confidence value for each amplitude component
cTol = 0.01;
pCN = 0.95;
ppfVal = ones(rows(FINALAMP),3);
%for count = 1:rows(ppfVal)
%	for count2 = 1:3
%		ppfVal(count,count2) = riceppf(FINALAMP(count,count2+1),FINALERR(count,count2+1),pCN,cTol,FINALAMP(count,count2+1) + 2*FINALERR(count,count2+1));
%	endfor
%endfor


%Thermal noise limit calculation for torque and g_{B-L}
thNoise = thermalNoise(FINALAMP(:,1),kappa,Q,Temp,f0,aCN,rows(checkLength),isExternal);
gLim = torqueToGBL(thNoise);

if(saveData)
	fileName = ['run0160',ctime(time)([5:7,9:10,21:end-1]),'(',num2str(startFreq),'-',num2str(stopFreq),')'];
	save(fileName,'FINALAMP','FINALERR','gLim');
	disp(['Saved run to file ',fileName]);
endif

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
%! disp('transferFunction');
%! test transferFunction
%! disp('twoPointTransfer');
%! test twoPointTransfer
%! disp('ampToPower');
%! test ampToPower
%! disp('preCalcComponents');
%! test preCalcComponents
%! disp('testerFile');
%! fflush(stdout);
