t = 1:9500; t=t';
Amp = 1;
f = 2*pi*(abs(randn)/10)
fData = [t,Amp.*sin(f.*t)];
%mean(fData)
%dX = [ones(rows(t),1),t];
%[b,s,r,err,cov] = ols2(fData(:,2),dX);
%fData(:,2) = fData(:,2) - dX*b;
%mean(fData)
%fflush(stdout);
dD = cell(1,1);
dD{1,1} = fData;
freqArray = ((0:(rows(t)/2))')./(rows(t));
[compAvg,compOut,errOut] = dispAmpTF(dD,freqArray,0,0,1,0);
fftOut = (1/rows(t)).*fft(fData(:,2));
fftOut = fftOut(1:(rows(fftOut)/2 + 1),:);

figure(1);
loglog(freqArray,abs(abs(fftOut)-abs(compAvg(:,1)/2)));
figure(2);
semilogx(freqArray,angle(fftOut)-angle(compAvg(:,1)));
%fakeSignal =(1/(1.5e6)).*ifft(-[compAvg(:,1);flip(conj(compAvg(2:end-1,1)))]);
f
%plot(t,[fData(:,2),fakeSignal]);
%legend('OG','NOT OG');
%for count = 1:rows(t)
%	if(abs(real(fakeSignal)(count) .- fData(count,2)) > 1)
%		[count,real(fakeSignal)(count) .- fData(count,2)]
%		fflush(stdout);
%	endif
%endfor
%assert(abs(real(fakeSignal) .- fData(:,2))< (Amp / 2));


