function [theta,tau,pointsRemove] = removeEarthquakes(data,torque,threshold,areaRemove)

%Sets up for removal of high torque points.
noTorqueDisplace = data;
noEarthTorque = torque;
count = 0;
%Goes through each point to check if torque is above threshold, sets points
%over threshold to zero
i=rows(torque);
while (i>0)
  if (abs(torque(i,2))>(threshold))
    %Removes all points within areaRemove of the earthquake point
    back = (i-areaRemove);
    forward = (i+areaRemove);
    if ((i-areaRemove)<0)
      back = 1;
    endif
    if ((i+areaRemove)>length(noTorqueDisplace))
      forward = length(noTorqueDisplace);
    endif
    noTorqueDisplace(back:forward,:) = [];
    noEarthTorque(back:forward,2) = 0;
    i = back;
  endif
  i = i - 1;
endwhile


%Returns edited array
pointsRemove = count;
theta = noTorqueDisplace;
tau = noEarthTorque;
endfunction