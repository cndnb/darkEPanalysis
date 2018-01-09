%Pendulum and balance parameters, in SI units:
I = 378/(1e7);                                                                    
f0 = 1.9338e-3;                                                                 
Q = 500000;                                                                     
T = 273+24;  
kappa = (2*pi*f0)^2 * I;



%Response function:
tau = ones(length(searchFreq),1);
for count = 1:length(searchFreq)
tau(count,1) = (1/kappa)/((1-(searchFreq(count)/f0)^2)+((i*searchFreq(count))/(Q*f0)));
endfor

out = ones(length(searchFreq),3);
for count = 1:length(searchFreq)
  [BETA,COV] = powerFinder(O,searchFreq(count),50);
  out(count,1) = searchFreq(count);
  out(count,2) = abs(real(sqrt(BETA*BETA')/tau(count,1)));
  out(count,3) = abs(sqrt(COV*COV')/tau(count,1));
endfor
