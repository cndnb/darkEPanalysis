function [theta,tau] = removeEarthquakes(data,torque,threshold,areaRemove,showOut)
if (nargin != 5)
  error('[theta,tau] = removeEarthquakes(data,torque,threshold,areaRemove,showOut)');
endif

noEarthTorque = torque;
noEarthquakes = data;

subThresh = torque(:,2) .- threshold.*ones(rows(torque),1);
onlyPositive = subThresh./(-abs(subThresh)) - 1;
eqInd = find(onlyPositive);

for count = 1:rows(eqInd)	
	if (showOut)
		count
		fflush(stdout);
	endif
	back = eqInd(count) - areaRemove(1,1);
	forward = eqInd(count) + areaRemove(1,2);
	if (back < 1)
		back = 1;
	endif
	if (forward > rows(data))
		forward = rows(data);
	endif
	noEarthquakes(back:forward,2) = Inf;
	noEarthTorque(back:forward,2) = 0;
endfor

%Removes points that were set to Inf in time series
removeInd = find(isinf(noEarthquakes(:,2)));
noEarthquakes(removeInd,:) = [];
		
%Returns edited arrays
theta = noEarthquakes;
tau = noEarthTorque;

endfunction

%!test
%! t = 1:10000; t=t';
%! fData = [t,sin((2*pi*1/100).*t)];
%! torqueData = [t,zeros(rows(t),1)];
%! areaRemove = 500;
%! threshold = .99;
%! [rTime,rTorque] = removeEarthquakes(fData,torqueData,threshold,areaRemove,0);
%! assert(rTime,fData);
%! pointEarthquake = 5001;
%! torqueData(pointEarthquake,2) = 1;
%! [rTime,rTorque] = removeEarthquakes(fData,torqueData,threshold,areaRemove,0);
%! assert(rTorque,[t,zeros(rows(torqueData),1)]);
%! cTime = fData; cTime(pointEarthquake-areaRemove:pointEarthquake+areaRemove,:) = [];
%! assert(rTime,cTime);
