freq = [1e-4,1e-3,1e-2]';
driftFix = data;
figure(9);
hold on;
for count = 1:rows(freq)
plot(driftFix(:,1),frequencyVariance(driftFix,freq(count),50));
endfor
hold off;
legend(num2str(freq(1)),num2str(freq(2)),num2str(freq(3)));%,num2str(freq(4)));

figure(10);
plot(driftFix(:,1),driftFix(:,2))