function [theta,tau] = removeEarthquakes(data,torque,threshold,areaRemove)
%This method works as long as areaRemove is less than one day.


dayLength = 86400; %seconds
numDays = floor(rows(data)/dayLength);
dataDivisions = cell(numDays,1);

for count = 0:numDays-1
  dataDivisions(count+1,1) = data((count*dayLength+1):((count+1)*dayLength),:);
endfor

%Prepares torque array to have points removed.
noEarthTorque = torque;

%Goes through each point to check if torque is above threshold, sets points
%over threshold to zero
for count = 1:numDays
  aMatrix = dataDivisions{count,1};
  indS = aMatrix(rows(aMatrix),1);
  while (indS>aMatrix(1,1)-1)
    if (abs(torque(indS,2))>(threshold)) %If torque at time exceeds threshold
      %Removes all points within areaRemove of the earthquake point
      
      %Displacement points removal
      back = ((indS-aMatrix(1,1))-areaRemove);
      forward = ((indS-aMatrix(1,1))+areaRemove);
      if (back < 1) %If removal area leaks into previous day
        try %Attempt to remove points from previous day
          subMatrix = dataDivisions{count - 1,1}; %Accesses next day
          subMatrix(rows(subMatrix)+back:rows(subMatrix),:) = []; %Back would be negative, so index = rows(subMatrix)+back
          dataDivisions{count - 1,1} = subMatrix;
          back = 1;
        catch %Otherwise matrix is the first day, do nothing.
          back = 1;
        end_try_catch
      endif
      if (forward > (rows(aMatrix))) %If removal area leaks into next day
        try
          subMatrix = dataDivisions{count + 1,1};
          subMatrix(1:(forward - rows(aMatrix)),:) = [];
          dataDivisions{count + 1,1} = subMatrix;
          forward = rows(aMatrix);
        catch
          forward = rows(aMatrix);
        end_try_catch
      endif
      aMatrix(back:forward,:) = [];
      dataDivisions{count,1} = aMatrix;
      
      %Removal on torque array
      back = indS - areaRemove;
      forward = indS + areaRemove;
      if (back < 1)
        back = 1;
      endif
      if (forward > rows(noEarthTorque))
        forward = rows(noEarthTorque);
      endif 
      noEarthTorque(back:forward,2) = 0;
      indS = back;
    endif
    indS = indS - 1;
  endwhile
endfor

%Returns edited arrays
theta = dataDivisions;
tau = noEarthTorque;

endfunction