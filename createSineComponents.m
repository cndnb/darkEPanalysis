function ret = createSineComponents(t,f)
  omegaSearch = 2*pi*f;
  omegaEarth = 2*pi*(1/86164.0916);
  
  X = ones(length(t),6);
  X(:,1)= sin(omegaSearch*t).*sin(omegaEarth*t);
  X(:,2)= sin(omegaSearch*t).*cos(omegaEarth*t);
  X(:,3)= cos(omegaSearch*t).*sin(omegaEarth*t);
  X(:,4)= cos(omegaSearch*t).*cos(omegaEarth*t);
  X(:,5)= sin(omegaSearch*t);
  X(:,6)= cos(omegaSearch*t);
  
  ret = X;
endfunction