function [theta,tau] = removeEarthquakes(data,torque,threshold,areaRemove,daysInclude)
%This method works as long as areaRemove is less than one day.
if (!exist('testing'))
  testing = 0;
endif

dayLength = 86164; %seconds
maxDays = ceil(data(rows(data),1)/dayLength);
if (daysInclude == 0) %Fits the maximum number of days in the data
  numDays = maxDays;
elseif (daysInclude > maxDays)
  numDays = maxDays;
else
  numDays = daysInclude;
endif

%Initialize accumulation array
dataDivisions = cell(maxDays,1);

endVal = data(end,1);
dayCount = 1;


modCount = mod(data(:,1),86400);
previousDay = 1;
lastNum = 1;
for secCount = 1:rows(data)
  secCount
  fflush(stdout);
  if(modCount(secCount) < lastNum)
    dataDivisions{dayCount,1} = data(previousDay:secCount - 1,:);
    previousDay = secCount;
    dayCount = dayCount + 1;
  endif
  lastNum = modCount(secCount);
endfor
dataDivisions{dayCount+1,1} = data(previousDay:rows(data),:);

%secCount = 1;
%while (dayCount <= numDays)
%	try
%		if(data(secCount,1) <= dayLength*dayCount)
%			dataDivisions{dayCount,1} = [dataDivisions{dayCount,1};data(secCount,:)];
%			secCount = secCount + 1;
%		else 
%			dayCount = dayCount + 1;
%		endif
%	catch
%		dayCount = dayCount + 1;
%	end_try_catch
%endwhile



%Prepares torque array to have points removed.
noEarthTorque = torque;

%Goes through each point to check if torque is above threshold, sets points
%over threshold to zero
count = maxDays;
while(count > 0)
	aMatrix = dataDivisions{count,1};
	if(rows(aMatrix) == 0)
		dataDivisions(count,:) = [];
	else
		indS = aMatrix(rows(aMatrix),1);
		noVal = 0;
		while (indS > aMatrix(1,1))
			if (abs(torque(indS,2))>(threshold)) %If torque at time exceeds threshold    
				%Removes all points within areaRemove of the earthquake point
				if (testing)
					count
					indS
					fflush(stdout);
				endif
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
  	endif
	count = count - 1;
endwhile

if(numDays > rows(dataDivisions))
  numDays = rows(dataDivisions);
endif

%Returns edited arrays
theta = dataDivisions(1:numDays,:);
tau = noEarthTorque;

endfunction
