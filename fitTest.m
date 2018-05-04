designX = createSineComponents(driftFix(:,1),1e-3);
[sChkBeta, sChkSigma, sChkR, sChkErr, sChkCov] = ols2(driftFix(:,2),designX);
plot(driftFix(:,1),driftFix(:,2),driftFix(:,1),designX*sChkBeta);
legend('data','fit');
ylabel('displacement (radians)')
xlabel('time (s)')