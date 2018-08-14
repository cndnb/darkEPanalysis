% Bandstops a specific frequency using a two-point filter and interpolation.
%
% The output is truncated in time by one period. 
% The timing is symmetric. f(t) = ( f(t-T/2) + f(t+T/2) )/2.0
% 
% Timing preserves the original sample cadence and phase.

function outData = torsionFilter(time, data, period)

%	outData = [];

%	for ctr = 1:rows(time)

%		t = time(ctr);
%		symmetricTimes = [ t - period/4.0 ; t + period/4.0];

%		if ( symmetricTimes(1)  > time( 1 ) & ...
%		     symmetricTimes(2)  < time(end) )

%			y2 = interp1( time, data, symmetricTimes, "linear");

%			outData = [outData ; t , mean(y2) ];
%		end
%	end

	symmetricTimes = [time - period/4.0, time + period/4.0];

	y1 = interp1(time, data, symmetricTimes(:,1),"linear");
	y2 = interp1(time, data, symmetricTimes(:,2),"linear");


	outData = [time mean( [y1 y2] ,2) ];

	outData = outData( !isnan(outData(:,2)) , :);
end

%!test
%!
%! t = 1:10000;
%! t = t';
%!
%! s = sin(2*pi*t/1000.3324);
%! o = torsionFilter(t,s, 1000.3324);
%! assert( max( abs( o(:,2) ) < 1e-8));
