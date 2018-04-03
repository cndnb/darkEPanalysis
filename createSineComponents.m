function ret = createSineComponents(t,f)
  %Simple statement of usage
  if (nargin != 2)
    usage ("X = createSineComponents(t,f)");
  endif
  
  omegaSearch = 2*pi*f;
  omegaEarth = 2*pi*(1/86164.0916);
  
  X = ones(length(t),6);
  %Perpendicular to X
  X(:,1)= sin(omegaSearch.*t).*sin(omegaEarth.*t);
  X(:,3)= cos(omegaSearch.*t).*sin(omegaEarth.*t);
  
  %Parallel to X
  X(:,2)= sin(omegaSearch.*t).*cos(omegaEarth.*t);
  X(:,4)= cos(omegaSearch.*t).*cos(omegaEarth.*t);
  
  %Z component
  X(:,5)= sin(omegaSearch.*t);
  X(:,6)= cos(omegaSearch.*t);
  
  ret = X;
endfunction

%!test
%! G = createSineComponents(0,1);
%! assert(G  == [ 0,0,0,1,0,1] );

%!test
%! b=[1,2,3,4,5,6,7,8,9,10];
%! b=b';
%! b = 10000.*b;
%! X = createSineComponents(b,pi);
%! assert (rank(X'*X) == 6)