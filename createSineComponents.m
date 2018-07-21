function [dZ,dPeX,dPaX] = createSineComponents(timeData,f)
  %Simple statement of usage
  if (nargin != 2)
    usage ("X = createSineComponents(t,f)");
  endif
  

  f0 = 1.9338e-3;
  omegaSearch = 2*pi.*f;
  omegaEarth = 2*pi*(1/86164.0916);
  
  %Creates 3-D array, 3rd dimension is the search frequency. Dim 1 and 2
  %is the design matrix for each frequency.
  dZ = ones(length(timeData),2);
  dPeX = ones(length(timeData),2);
  dPaX = ones(length(timeData),2);
  
  %Perpendicular to X
  dPeX(:,1)= sin(omegaSearch.*timeData).*sin(omegaEarth.*timeData);
  dPeX(:,2)= cos(omegaSearch.*timeData).*sin(omegaEarth.*timeData);
  
  %Parallel to X
  dPaX(:,1)= sin(omegaSearch.*timeData).*cos(omegaEarth.*timeData);
  dPaX(:,2)= cos(omegaSearch.*timeData).*cos(omegaEarth.*timeData);
  
  %Z component
  dZ(:,1)= sin(omegaSearch.*timeData);
  dZ(:,2)= cos(omegaSearch.*timeData);
  
  %Resonant frequency component
  %X(:,7) = sin((2*pi*f0).*timeData);
  %X(:,8) = cos((2*pi*f0).*timeData);
  
  %Daily mondulation component
  %X(:,9) = sin(omegaEarth.*timeData);
  %X(:,10) = cos(omegaEarth.*timeData);
  
  %Drift component
  %X(:,9) = timeData;
  
  %Constant offset component
  %X(:,3) = ones(rows(timeData),1);

endfunction

%!test
%! G = createSineComponents(0,1);
%! assert(G  == [ 0,0,0,1,0,1,0,1,0,1] );

%!test
%! b=1:10000;
%! b=b';
%! X = createSineComponents(b,pi);
%! assert (rank(X'*X) == columns(X))

%!test
%! b=1:10000;
%! b=b';
%! X = createSineComponents(b,pi);
%! assert (b == X(:,(columns(X)-1)))
