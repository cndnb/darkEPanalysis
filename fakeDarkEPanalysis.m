function [AMPOUT,AMPVAR] = fakeDarkEPanalysis(data,chunkSize, jump, startFreq, stopFreq,threshold)
%%%%%%%%%%%%%%%%%%%% PROBLEM LAYOUT & CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Damped oscillator differential equation: x'' + (2wZ) x' + (w^2) x = 0
%Where w is the undamped angular frequency, and Z = 1/2Q. Oscillation is underdamped,
%therefore solution is of the form x = e^(-wZt)[A sin(Ct) + B cos(Ct)]
%C = sqrt[w^2-(wZ)^2] = w*sqrt[1-Z^2]

%torque = kappa*theta + rotI*(d^2 theta/dt^2)

%Pendulum and balance parameters, in SI units:
I = 378/(1e7);                                                                    
f0 = 1.9338e-3;                                                                 
Q = 500000;                                                                     
T = 273+24;  
kappa = (2*pi*f0)^2 * I;

%Variables important to fitting
endCount = (stopFreq-startFreq)/(jump*(1/length(data)));
%Creates plotting array
ampFreq = ones((stopFreq-startFreq)/(jump*(1/length(data))),7);
%Creates error array
ampVar = ones((stopFreq-startFreq)/(jump*(1/length(data))),7);

for count = 1:endCount
[BETA, COV] = powerFinder(data,...
(startFreq+(count*jump*(1/length(data)))),chunkSize,threshold);


ampFreq(count,:) = [(startFreq+(count*jump*(1/length(data)))),BETA];
ampVar(count,:) = [(startFreq+(count*jump*(1/length(data)))),COV];

endfor
AMPOUT = ampFreq;
AMPVAR = ampVar;
endfunction