%%%%%%%%%%%%%%%%%%%% UNIT TESTS/PRELIMINARY CHECKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testing = 0;

if (testing)
  'dispAmpTF'
  test dispAmpTF
  'createSineComponents'
  test createSineComponents
  'frequencyVariance'
  test frequencyVariance
  'specFreqAmp'
  test specFreqAmp
  'weightedOLS'
  test weightedOLS
  'transferFunction'
  test transferFunction
  'ampToPower'
  test ampToPower

  fflush(stdout);
  pause()
endif

%%%%%%%%%%%%%%%%%%%% PROBLEM LAYOUT & CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Damped oscillator differential equation: x'' + (2wZ) x' + (w^2) x = 0
%Where w is the undamped angular frequency, and Z = 1/2Q. Oscillation is underdamped,
%therefore solution is of the form x = e^(-wZt)[A sin(Ct) + B cos(Ct)]
%C = sqrt[w^2-(wZ)^2] = w*sqrt[1-Z^2]

%torque = kappa*theta + rotI*(d^2 theta/dt^2)

%Pendulum and balance parameters, in SI units:
global I = 378/(1e7);                                                                    
global f0 = 1.9338e-3;                                                                 
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
linearColumn = 0;%numBETAVal - 1;

%Multiples of smallest usable frequency between amplitude points
jump = 1;
%Start of frequency scan
startFreq = 1e-3;
%End frequency scan
stopFreq = 1e-2;

%If weighted
%1 for weightedOLS, 0 for ols2
fitIsWeighted = 0;

%Number of days in the data considered
daysInclude = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%How many periods of the specific frequency are included in weighted error fit
chunkSize = 10;

pkg load signal;

if (!exist('d'))
  importfakeDarkEP
  %f = fir1(10000,0.01,'high');
  %F = filter(f,1,d(:,2));
  %weightVal = abs(F - mean(F)) + 1e-9;
  weightVal = resonanceVariance(d,chunkSize);
  d(:,3) = weightVal;
  newD = [d(195000:235000,:);d(370000:410000,:)];
endif
omegaEarth = 2*pi*(1/86164.0916);
t = (1:rows(newD))';
X = [ones(rows(t),1)];%,t,sin(omegaEarth.*t),cos(omegaEarth.*t)];
[DFB,DFS,DFR,DFERR,DFCOV] = ols2(newD(:,2),X);
preDF = [newD(:,1),newD(:,2) - X*DFB];%,newD(:,3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%% EARTHQUAKE REMOVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculates the torque at each point, puts into an array for analysis
calcTorque = torque(newD, I, kappa);

%This is the value above which torque is considered an earthquake
threshold = 1e-13 + mean(calcTorque(3:(1e6-2),2));
%Number of seconds around a large torque that will be removed
areaRemove = 10000;


%returns torques set to zero at earthquakes in a matrix, 
%driftFix = data divided into days and earthquake points removed
%driftFix{day,1} = [seconds, displacement amplitude]
%Full length is length of the data in seconds from start to stop, before
%earthquake removal
[driftFix,editTorque] = removeEarthquakes(preDF,calcTorque,threshold,areaRemove,daysInclude);

fullLength = rows(driftFix{1,1});

if (testing)
  %Makes plotting more simple
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
  
  figure(1);
  plot(d(:,1),d(:,2));
  title('Original Data');
  xlabel('Time (s)');
  ylabel('Displacement (rad)');
  
  figure(2);
  plot(fullData(:,1),fullData(:,2));
  title('Earthquakes Removed');
  xlabel('Time (s)');
  ylabel('Displacement (rad)');
  
  %Checking that FFT of torque has no peaks
  check = psd(editTorque(:,1),editTorque(:,2)-mean(editTorque(:,2)));
  figure(5);
  loglog(check(:,1),check(:,2));
  title('Torque FFT without earthquakes');
  xlabel('Time (s)');
  ylabel('Torque (N m)');
  
  %Checking the threshold level
  figure(6);
  plot(calcTorque(3:(1e6 - 2),1),calcTorque(3:(1e6 - 2),2),calcTorque(3:(1e6 - 2),1),...
  threshold.*ones(length(calcTorque)-4,1));
  title('Threshold plotted on torque');
  xlabel('Time (s)');
  ylabel('Torque (N m)');
endif

%%%%%%%%%%%% AMPLITUDE(TIME) => AMPLITUDE(FREQUENCY) CONVERSION  %%%%%%%%%%%%%%%%

%Fits design matrix at each frequency over a given number of data bins
%ampFreq is the average of each amplitude at each frequency over all bins
%ampError is the standard deviation of each amplitude at each frequency over all bins

%endCount is the #rows of frequency matrix
%endCout = total frequency band divided by the smallest frequency jump
%Integer so that it can be used for indexing
endCount = floor((stopFreq-startFreq)/(jump*(1/fullLength)))+1
fflush(stdout);

freqArray = ones(endCount,1);
  
%Assigns frequency values for the first column of the frequency and error arrays
for count = 1:endCount
  freqArray(count,1) = (startFreq+((count-1)*jump*(1/fullLength))); %fullLength passed before earthquakes removed
endfor
  
%[ampFreq,ampError] = dispAmpTF(driftFix,freqArray,endCount,linearColumn,fitIsWeighted,1);
[preAvgZ,preAvgPerpX,preAvgParaX] = dispAmpTF(driftFix,freqArray,endCount,linearColumn,fitIsWeighted,1);

%Correction for fake dataset
freq = 9e-3;
phaseCorrection = ones(rows(freqArray),1,daysInclude);
for count = 1:rows(freqArray)
	for dayCount = 1:daysInclude
		[dZ,dPeX,dPaX] = createSineComponents(driftFix{dayCount,1}(:,1),freqArray(count));
		[b1,s1,r1,err1,cov1] = ols2(driftFix{dayCount,1}(:,2),dZ);
		phaseCorrection(count,1,dayCount) = angle(b1(2,1) + i.*b1(1,1));
	endfor
endfor
compOut = preAvgZ(:,2,:) + i.*preAvgZ(:,1,:);
compOut = compOut.*phaseCorrection;
assert(angle(compOut),zeros(size(compOut)),2*eps);
preAvgZ = [real(compOut),imag(compOut)];

%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSION TO TORQUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Creates array for final data
FINALAMP = ones(rows(freqArray),3);%,5);
FINALERR = ones(rows(freqArray),3);%,4);

%Sums in quadrature amplitudes to find single value for each coordinate direction,
%Divides by the transfer function to find the torque amplitude for each frequency
[FINALAMP, FINALERR] = ampToPower(preAvgZ,preAvgPerpX,preAvgParaX,freqArray,kappa,f0,Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plots torque power as a function of frequency
figure(1);
loglog(FINALAMP(:,1),FINALAMP(:,2));
%hold on;
%loglog([FINALAMP(:,1),FINALAMP(:,1)],[FINALAMP(:,2) - FINALERR(:,2),FINALAMP(:,2) + FINALERR(:,2)],'-r');
%hold off;
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency parallel to Z');

figure(2);
loglog(FINALAMP(:,1),FINALAMP(:,3));
%hold on;
%loglog([FINALAMP(:,1),FINALAMP(:,1)],[FINALAMP(:,3) - FINALERR(:,3),FINALAMP(:,3) + FINALERR(:,3)],'-r');
%hold off;
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency perpendicular to gamma');

figure(3);
loglog(FINALAMP(:,1),FINALAMP(:,4));
%hold on;
%loglog([FINALAMP(:,1),FINALAMP(:,1)],[FINALAMP(:,4) - FINALERR(:,4),FINALAMP(:,4) + FINALERR(:,4)],'-r');
%hold off;
xlabel('Frequency (Hz)');
ylabel('Torque (N m)');
title('Torque vs frequency parallel to gamma');

%figure(4);
%loglog(FINALAMP(:,1),FINALAMP(:,2),FINALAMP(:,1),FINALAMP(:,3),FINALAMP(:,1),FINALAMP(:,4),FINALAMP(:,1),FINALAMP(:,5));
%legend('Parallel to gamma','Perpendicular to gamma','Z component','Sum signal');
%xlabel('Frequency (Hz)');
%ylabel('Torque (N m)');
%title('Torque vs frequency');
