function [theta,tau] = removeEarthquakes(data,torque,externalRemove,areaRemove,showOut)
if (nargin != 5)
  error('[theta,tau] = removeEarthquakes(data,torque,externalRemove,areaRemove,showOut)');
endif

noEarthTorque = torque{1,1};
noEarthquakes = data;

subThresh = torque{1,1}(:,2) .- torque{1,2}.*ones(rows(torque{1,1}),1);
onlyPositive = subThresh./(-abs(subThresh)) - 1;

%Indexes to be removed from earthquakes
indRemT = find(onlyPositive);
%Indexes to be removed from external triggers
indRemX = find(externalRemove);

%All indexes to be removed in order
eqInd = sort([indRemT;indRemX]);

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
	noEarthquakes(back:forward,2) = NaN;
	noEarthTorque(back:forward,2) = 0;
endfor


%Removes points that were set to Inf in time series
removeInd = find(isnan(noEarthquakes(:,2)));
noEarthquakes(removeInd,:) = [];
		
%Returns edited arrays
theta = noEarthquakes;
tau = noEarthTorque;

endfunction

%!test
%! t = 1:10000; t=t';
%! fData = [t,sin((2*pi*1/100).*t)];
%! torqueData = cell(1,2);
%! torqueData{1,1} = [t,zeros(rows(t),1)];
%! torqueData{1,2} = .99; %Threshold
%! areaRemove = [500,500];
%! externalRemove = zeros(rows(fData),1);
%! [rTime,rTorque] = removeEarthquakes(fData,torqueData,externalRemove,areaRemove,0);
%! assert(rTime,fData);
%! pointEarthquake = 5001;
%! torqueData{1,1}(pointEarthquake,2) = 1;
%! [rTime,rTorque] = removeEarthquakes(fData,torqueData,externalRemove,areaRemove,0);
%! assert(rTorque,[t,zeros(rows(torqueData{1,1}),1)]);
%! cTime = fData; cTime(pointEarthquake-areaRemove(1,1):pointEarthquake+areaRemove(1,2),:) = [];
%! assert(rTime,cTime);
