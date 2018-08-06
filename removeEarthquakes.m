function [theta,tau] = removeEarthquakes(data,torque,threshold,areaRemove)
if (nargin != 4)
	error('[theta,tau] = removeEarthquakes(data,torque,threshold,areaRemove)');
endif

%Prepares torque and data arrays to have points removed.
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
		noEarthquakes(back:forward,2) = 0;
		noEarthTorque(back:forward,2) = 0;
	catch
	end_try_catch
	count = count - 1;
endwhile

removeInd = find(noEarthquakes(:,2)./noEarthquakes(:,2) - 1);

noEarthquakes(removeInd,:) = [];
		
%Returns edited arrays
theta = noEarthquakes;
tau = noEarthTorque;

endfunction
