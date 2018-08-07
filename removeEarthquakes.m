function [theta,tau] = removeEarthquakes(dataDivisions,torque,threshold,areaRemove,daysInclude)
%This method works as long as areaRemove is less than one day.
if (!exist('testing'))
  testing = 0;
endif

%Prepares torque array to have points removed.
noEarthTorque = torque;

%Goes through each point to check if torque is above threshold, sets points
%over threshold to zero
count = numDays;
while(count > 0)
	aMatrix = dataDivisions{count,1};
	if(rows(aMatrix) == 0)
		dataDivisions(count,:) = [];
	else
		indS = rows(aMatrix);
		noVal = 0;
		while (indS > 0)
			if (abs(torque(aMatrix(indS,1),2))>(threshold)) %If torque at time exceeds threshold    
				%Removes all points within areaRemove of the earthquake point
				if (testing)
					count
					indS
					fflush(stdout);
				endif
				%Displacement points removal
				back = (indS-areaRemove);
				forward = (indS+areaRemove);
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
