function [AMPOUT,AMPVAR] = fakeDarkEPanalysis(data,chunkSize, jump, startFreq, endCount)
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
%Creates plotting array
ampFreq = ones(endCount,7);
%Creates error array
ampVar = ones(endCount,7);
for count = 1:endCount
  count
  fflush(stdout);
[BETA,COV] = specFreqPower(data,...
(startFreq+((count-1)*jump*(1/rows(data)))),chunkSize);
ampFreq(count,:) = [(startFreq+((count-1)*jump*(1/rows(data)))),BETA];
ampVar(count,:) = [(startFreq+((count-1)*jump*(1/rows(data)))),COV];

endfor
AMPOUT = ampFreq;
AMPVAR = ampVar;
endfunction