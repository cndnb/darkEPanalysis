function answer=psdw(timebase,data)

	if(nargin==1)
		data=timebase(:,2);
		timebase=timebase(:,1);
	end

	if(nargin>2)
		usage("No more than two arguments allowed!");
	end
	
	timebase=timebase-timebase(1);
	sampleTime=timebase(2)

	answer=timebase;

       %Windowing - Welch
 
       N = length(data);
 
       window = (1:N)';
 
       window = 1 - ( 2* (window - N/2)./N).^2;
 
       wdata = window.*data;

 
	df=fft(wdata);
	answer(:,2)=df.*transpose(df')/length(data)^2;
	%answer(:,1)=1./(timebase(length(timebase)/2)+1-timebase);
	answer(:,1)=timebase.*1/(timebase(length(timebase))-timebase(1))/(timebase(2)-timebase(1));

	%Parseval's theorem that holds here: 

	%dataPower=sum(data.*data)*timebase(2)
	dataPower=sum(data.*data)/length(data);
%	psdPower=sum(answer(:,2))*answer(2,1)
	psdPower=sum(answer(:,2));

	corf=dataPower/psdPower
	answer(:,2)=answer(:,2)*corf;

%	psdPower=sum(answer(:,2))/length(data)^2;

	%Make PSD one-sided <-- is there a problem with odd length data?
	answer = answer(1:ceil(length(answer)/2),:);
	answer(:,2) = 2*answer(:,2);

%	for i=1:length(answer)/2
%		answer(i,2)=answer(i,2)+answer(length(answer)+1-i,2);
%	end

	%If there is a problem with odd length data, it gets chopped here, methinks.
%	answer=answer(1:(length(data)/2),:);

%	answer(:,2)=answer(:,2)*2; %Corrects for truncating - fix with reflecting sum;

	% this should be the same as data power
	answerSum= sum(answer(:,2))
	% may not be the same as data power for windowed data
	psdPower
	dataPower
	answerCorf=sum(answer(:,2))/psdPower

	%Normalization: At this stage, the mean squared amplitude of the data is equal to the sum of the psd.

	%Makes different length psds lie on top of one another. The sum of the psd is now equal to the sum of the squares of the data; May still need a factor of two or something.
	answer(:,2)=answer(:,2).*length(data)*sampleTime;






endfunction

