function rtn = createSineComponents(timeData,f,cM)
  if (nargin != 3)
    error("X = createSineComponents(t,f,seattleLat,seattleLong,compassDir,startTime)");
  endif
  if (size(cM) != [rows(timeData),3])
    error("constant matrix has incorrect size");
  endif
  

  f0 = 1.9338e-3;
  omegaSearch = 2*pi.*f;
  omegaEarth = 2*pi*(1/86164.0916);
  oED = 2*pi*(1/86400);

    
  X = ones(rows(timeData),12);  

  OMEGA = omegaSearch.*timeData;
  resOMEGA = 2*pi*f0.*timeData;
  eOMEGA = omegaEarth.*timeData;
  sCCol = [sin(OMEGA),cos(OMEGA)];  


  %Z component
  X(:,1:2) = sCCol.*cM(:,1); 

  %Perpendicular to X
  X(:,3:4) = sCCol.*cM(:,2);  

  %Parallel to X
  X(:,5:6) = sCCol.*cM(:,3);
  
  %Daily mondulation component
  X(:,7) = sin(eOMEGA);
  X(:,8) = cos(eOMEGA);
  
  %Drift component
  X(:,9) = timeData;
  
  %Constant offset component
  X(:,10) = ones(rows(timeData),1);
 
  %Resonant frequency component
  X(:,11) = sin(resOMEGA);
  X(:,12) = cos(resOMEGA);

  rtn = X;
endfunction

%!test
%! cM = preCalcComponents(0,0,0,0,0);
%! G = createSineComponents(0,1,cM);
%! assert(G  == [ 0,1,0,0,0,0,0,1,0,1,0,1] );

%!test
%! b=1:86400;
%! b=b';
%! f = 1.933e-3;
%! cM = preCalcComponents(b,pi/4,pi/4,0,0);
%! X = createSineComponents(b,f,cM);
%! assert (rank(X'*X) == columns(X))

