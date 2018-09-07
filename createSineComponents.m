function rtn = createSineComponents(timeData,f)
  if (nargin != 2)
    usage ("X = createSineComponents(t,f)");
  endif
  

  f0 = 1.9338e-3;
  omegaSearch = 2*pi.*f;
  omegaEarth = 2*pi*(1/86164.0916);
  oED = 2*pi*(1/86400);
  
  %Creates 3-D array, 3rd dimension is the search frequency. Dim 1 and 2
  %is the design matrix for each frequency.
  %dZ = ones(length(timeData),2);
  %dPeX = ones(length(timeData),2);
  %dPaX = ones(length(timeData),2);
  X = ones(rows(timeData),6);  
  
  %Z component
  X(:,1)= sin(omegaSearch.*timeData);
  X(:,2)= cos(omegaSearch.*timeData);
  %dZ(:,1) = sin(omegaSearch.*timeData);
  %dZ(:,2) = cos(omegaSearch.*timeData);

  %Perpendicular to X
  X(:,3) = sin(omegaSearch.*timeData).*sin(omegaEarth.*timeData);
  X(:,4) = cos(omegaSearch.*timeData).*sin(omegaEarth.*timeData);
  %dPeX(:,1) = sin(omegaSearch.*timeData).*sin(omegaEarth.*timeData);
  %dPeX(:,2) = cos(omegaSearch.*timeData).*sin(omegaEarth.*timeData);
  
  %Parallel to X
  X(:,5) = sin(omegaSearch.*timeData).*cos(omegaEarth.*timeData);
  X(:,6) = cos(omegaSearch.*timeData).*cos(omegaEarth.*timeData);
  %dPaX(:,1) = sin(omegaSearch.*timeData).*cos(omegaEarth.*timeData);
  %dPaX(:,2) = cos(omegaSearch.*timeData).*cos(omegaEarth.*timeData);
  
  %Daily mondulation component
  X(:,7) = sin(omegaEarth.*timeData);
  X(:,8) = cos(omegaEarth.*timeData);
  
  %Drift component
  X(:,9) = timeData;
  
  %Constant offset component
  X(:,10) = ones(rows(timeData),1);
  
  %Resonant frequency component
  X(:,11) = sin((2*pi*f0).*timeData);
  X(:,12) = cos((2*pi*f0).*timeData);

  rtn = X;
endfunction

%!test
%! G = createSineComponents(0,1);
%! assert(G  == [ 0,1,0,0,0,1] );

%!test
%! b=1:10000;
%! b=b';
%! X = createSineComponents(b,pi);
%! assert (rank(X'*X) == columns(X))

%!test
%! b=10000:20000;
%! b=b';
%! X = createSineComponents(b,pi);
%! f = 2*pi*pi;
%! fE = 2*pi*(1/86164.0916); 
%! checkX = [sin(f.*b),cos(f.*b),sin(f.*b).*sin(fE.*b),cos(f.*b).*sin(fE.*b),sin(f.*b).*cos(fE.*b),cos(f.*b).*cos(fE.*b)];
%! assert(X,checkX);
