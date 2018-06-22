'importing data'
fflush(stdout);

clear driftFix;
clear fullLength;

if (!exist('testing'))
  testing = 0;
endif

%imports data
d = load('fakeDarkEPAugust92017.dat');

%%%%%%%%%%%%%%%%%%%%%%%% EARTHQUAKE REMOVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculates the torque at each point, puts into an array for analysis
calcTorque = torque(d, I, kappa);

%This is the value above which torque is considered an earthquake
threshold = 1e-13 + mean(calcTorque(3:(1e6-2),2));
%Number of seconds around a large torque that will be removed
areaRemove = 10000;
%Number of days in the data considered
daysInclude = 3;

%returns torques set to zero at earthquakes in a matrix, 
%driftFix = data divided into days and earthquake points removed
%driftFix{day,1} = [seconds, displacement amplitude]
%Full length is length of the data in seconds from start to stop, before
%earthquake removal
[driftFix,editTorque,fullLength] = removeEarthquakes(d,calcTorque,threshold,areaRemove,daysInclude);






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

'done'
fflush(stdout);