if (!exist('d'))
	importfakeDarkEP
endif
I = 378/(1e7);
f0 = 1.9338e-3;
kappa = (2*pi*f0)^2 * I;

calcTorque = torque(d,I,kappa);

%This is the value above which torque is considered an earthquake
threshold = 1e-13 + mean(calcTorque(3:(end-2),2));
%Number of seconds around a large torque that will be removed
areaRemove = 10000;

%Returns array with torque points above threshold zeroed
[noEQ,noTorque] = removeEarthquakes(d,calcTorque,threshold,areaRemove);

figure(1);
check = psd(noTorque(:,1),noTorque(:,2)-mean(noTorque(:,2)));
loglog(check(:,1),check(:,2));

figure(2);
plot(d(:,1),d(:,2),noEQ(:,1),noEQ(:,2));
legend('OG','REMOVED');
