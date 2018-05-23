function ret = createSineComponents(t,f)
  %Simple statement of usage
  if (nargin != 2)
    usage ("X = createSineComponents(t,f)");
  endif
  
  f0 = 1.9338e-3;
  omegaSearch = 2*pi.*f;
  omegaEarth = 2*pi*(1/86164.0916);
  
  %Creates 3-D array, 3rd dimension is the search frequency. Dim 1 and 2
  %is the design matrix for each frequency.
  X = ones(length(t),12);
  
  %Perpendicular to X
  X(:,1)= sin(omegaSearch.*t).*sin(omegaEarth.*t);
  X(:,3)= cos(omegaSearch.*t).*sin(omegaEarth.*t);
  
  %Parallel to X
  X(:,2)= sin(omegaSearch.*t).*cos(omegaEarth.*t);
  X(:,4)= cos(omegaSearch.*t).*cos(omegaEarth.*t);
  
  %Z component
  X(:,5)= sin(omegaSearch.*t);
  X(:,6)= cos(omegaSearch.*t);
  
  %Resonant frequency component
  X(:,7) = sin((2*pi*f0).*t);
  X(:,8) = cos((2*pi*f0).*t);
  
  %Daily mondulation component
  X(:,9) = sin(omegaEarth.*t);
  X(:,10) = cos(omegaEarth.*t);
  
  %Drift component
  X(:,11) = t;
  
  %Constant offset component
  X(:,12) = ones(rows(t),1);
  
  ret = X;
endfunction

%!test
%! G = createSineComponents(0,1);
%! assert(G  == [ 0,0,0,1,0,1,0,1,0,1,0,1] );

%!test
%! b=1:10000;
%! b=b';
%! X = createSineComponents(b,pi);
%! assert (rank(X'*X) == columns(X))