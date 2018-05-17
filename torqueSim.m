function TimeSeries = torqueSim(time,  momInertia, kappa, Q, Temperature, ExtTor)
  if (nargin !=6)
    usage('TimeSeries = torqueSim(t,I,kappa,Q,T,tau)');
  endif
	%CODATA 2002
	k_B=1.3806505e-23;

	dt = mean(diff(time));

	%implements velocity damping only

	theta = 0;
	w = 0;

	TimeSeries = zeros(rows(time), 4);
  
	w0 = sqrt(kappa/momInertia);

	  TorN = sqrt(4*k_B * Temperature * kappa/Q / w0) * randn(rows(time),1);
  
	for ctr = 1:rows(time)
	
		%equation of motion:
		%w = w + (-kappa * theta - w * kappa / Q / w0 + TorN(ctr) + ExtTor(ctr)) * dt / momInertia;
    w = w + (-kappa * theta - w * kappa / Q / w0  + TorN(ctr) + ExtTor(ctr)) * dt / momInertia;
		%w = w + (-kappa * theta ) * dt / momInertia;
		theta = theta + w * dt;

		TimeSeries(ctr, :) = [time(ctr) theta w TorN(ctr)];
	end

end

%!test
%! darkEPTEST(1e5,100)
