function [theta,tau] = removeEarthquakes(data,torque,threshold,areaRemove)
if (nargin != 4)
  error('[theta,tau] = removeEarthquakes(data,torque,threshold,areaRemove)');
endif

noEarthTorque = torque;
noEarthquakes = data;

subThresh = torque(:,2) .- threshold.*ones(rows(torque),1);
onlyPositive = subThresh./(-abs(subThresh)) - 1;
eqInd = find(onlyPositive);

count = rows(eqInd);
while(count > 0)
	count
	fflush(stdout);
	back = eqInd(count) - areaRemove;
	forward = eqInd(count) + areaRemove;
	if (back < 1)
		back = 1;
	endif
	if (forward > rows(data))
		forward = rows(data);
	endif
	try
		noEarthquakes(back:forward,2) = Inf;
		noEarthTorque(back:forward,2) = 0;
	catch
	end_try_catch
	count = count - 1;
endwhile

%Removes points that were set to Inf in time series
removeInd = find(isInf(noEarthquakes(:,2)));
noEarthquakes(removeInd,:) = [];
		
%Returns edited arrays
theta = noEarthquakes;
tau = noEarthTorque;

endfunction
