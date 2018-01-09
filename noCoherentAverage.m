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


%Variables important to fitting
%How many periods of the specific frequency are included in error fit
chunkSize = 50;
%Multiples of smallest usable frequency between amplitude points
jump = 100;
%Start of frequency scan
startFreq = 1e-6;
%End frequency scan
stopFreq = 1e-1;
%How many points are included in the coherent average bins
dataCut = rows(driftFix);

%endCount is the #rows of frequency matrix
endCount = (stopFreq-startFreq)/(jump*(1/dataCut));
%Creates plotting array
ampFreq = zeros(endCount,2);
%Creates error array
ampError = zeros(endCount,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Assigns frequency values for the first column of the frequency and error arrays
for i = 1:endCount
  ampFreq(i,1) = (startFreq+((i-1)*jump*(1/dataCut)));
endfor
ampError(:,1) = ampFreq(:,1);


%Initializes the response function for each frequency of the amplitudes
tau = ones(rows(ampFreq),1);
for count = 1:rows(ampFreq)
tau(count,1) = (1/kappa)/((1-(ampFreq(count,1)/f0)^2)+((i*ampFreq(count,1))/(Q*f0)));
endfor
'start'
startFreq
'stop'
stopFreq

for i = 1:rows(ampFreq)
  ampFreq(i,1)
  fflush(stdout);
  [BETA,ERR]=powerFinder(driftFix,ampFreq(i,1),chunkSize);
  ampFreq(i,2) = abs(real(sqrt((BETA*BETA'))/tau(i)));
endfor

figure(1);
loglog(ampFreq(:,1),ampFreq(:,2));

%%Creates array for final data
%FINALAMP = ones(rows(ampFreq),2);
%FINALERR = ones(rows(ampError),2);
%for count = 1:rows(ampFreq)
%%Sums in quadrature the averaged sine/cosine amplitudes (inner product)
%FINALAMP(count,:) = [ampFreq(count,1),abs(real(sqrt(ampFreq(count,[3,5])*ampFreq(count,[3,5])')/tau(count))),...
%abs(real(sqrt(ampFreq(count,[2,4])*ampFreq(count,[2,4])')/tau(count))),...
%abs(real(sqrt(ampFreq(count,6:7)*ampFreq(count,6:7)')/tau(count))),...
%abs(real(sqrt(ampFreq(count,2:7)*ampFreq(count,2:7)')/tau(count)))];
%endfor

%%Plots torque power as a function of frequency
%figure(1);
%loglog(FINALAMP(:,1),FINALAMP(:,2),FINALAMP(:,1),FINALAMP(:,3),FINALAMP(:,1),FINALAMP(:,4),FINALAMP(:,1),FINALAMP(:,5));
%legend('Parallel to gamma','Perpendicular to gamma','z component','Sum');
%xlabel('Frequency (Hz)');
%ylabel('Torque (N*m)');
%title('Torque vs frequency');