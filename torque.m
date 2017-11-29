function ret = torque(data,I,kappa)
%Approximate derivative for use in the torque differential equation
secondDerivative = [data(:,1),[0;0;symmetricDerivative(symmetricDerivative(data(:,2)));0;0]];

%This is the differential equation for torque
ret = [data(:,1), kappa.*data(:,2) + I.*secondDerivative(:,2)];
endfunction