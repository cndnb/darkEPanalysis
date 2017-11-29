function [theta,tau] = removeEarthquakes(data,torque,threshold,areaRemove)

%Sets up for removal of high torque points.
noTorqueDisplace = data;
noEarthTorque = torque;
%Goes through each point to check if torque is above threshold, sets points
%over threshold to zero
for i=1:length(torque)
  if (abs(torque(i,2))>(threshold))
    %Removes all points within areaRemove of the earthquake point
    back = (i-areaRemove);
    forward = (i+areaRemove);
    if ((i-areaRemove)<0)
      back = 1;
    endif
    if ((i+areaRemove)>length(torque))
      forward = length(torque)
    endif
    noTorqueDisplace(back:forward,2) = 0;
    noEarthTorque(back:forward,2) = 0;
  endif
endfor
%Returns edited array
theta = noTorqueDisplace;
tau = noEarthTorque;
endfunction