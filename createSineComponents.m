function rtn = createSineComponents(timeData,f,cM,cS)
	if (nargin != 4)
		error("X = createSineComponents(t,f,orientationConstants,columnSelector)");
	endif
	if (size(cM) != [rows(timeData),3])
		error(strcat("createSineComponents: constant matrix has incorrect size:",num2str(size(cM))));
	endif 

	f0 = 1.9338e-3;
	omegaSearch = 2*pi.*f;
	omegaEarth = 2*pi*(1/86164.0916);
	oED = 2*pi*(1/86400);

    
	X = ones(rows(timeData),10);  

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
  
	%Drift component
	X(:,7) = timeData .- (timeData(1,1) - 1);
  
	%Constant offset component
	X(:,8) = ones(rows(timeData),1);
 
	%Resonant frequency component
	X(:,9) = sin(resOMEGA);
	X(:,10) = cos(resOMEGA);
	
	rtn = X(:,find(cS)); 
endfunction

%!test
%! cM = preCalcComponents(0,0,0,0,0);
%! G = createSineComponents(0,1,cM,ones(1,10));
%! assert(G  == [ 0,1,0,0,0,0,1,1,0,1]);

%!test
%! cM = preCalcComponents(0,0,0,0,0);
%! assert(createSineComponents(0,1,cM,[0,1,0,0,0,0,1,1,0,1]) == ones(1,4));

%!test
%! seattleLat = pi/4; seattleLong = pi/4; compassDir = pi/6; startTime = 0;f0 = 1.9338e-3;
%! t = (1:100)';
%! cM = preCalcComponents(t,seattleLat,seattleLong,compassDir,startTime);
%! f = rand(1)./10;
%! resOMEGA = 2*pi*f0.*t;
%! OMEGA = 2*pi*f.*t;
%! sCCol = [sin(OMEGA),cos(OMEGA)];
%! X = ones(100,10);
%! X(:,1:2) = sCCol.*cM(:,1); 
%! X(:,3:4) = sCCol.*cM(:,2);  
%! X(:,5:6) = sCCol.*cM(:,3);
%! X(:,7) = t .- (t(1,1) - 1);
%! X(:,8) = ones(rows(t),1);
%! X(:,9) = sin(resOMEGA);
%! X(:,10) = cos(resOMEGA);
%! 
%! for count = 1:10
%! 	columnSelector = zeros(1,10);
%! 	columnSelector(1,count) = 1;
%!	tX = createSineComponents(t,f,cM,columnSelector);
%!	assert(tX == X(:,count));
%! endfor

%!test
%! t = (1:100)';
%! cM = [ones(rows(t),1),zeros(rows(t),2)];
%! f = rand(1)./10;
%! columnSelector = [1,1,0,0,0,0,0,0,0,0];
%! X = createSineComponents(t,f,cM,columnSelector);
%! compX = [sin(2*pi*f.*t),cos(2*pi*f.*t)];
%! assert(X == compX);

%!test
%! b=1:86400;
%! b=b';
%! f = 1.933e-3;
%! f2 = 2*pi*(1/86164.0916);
%! cM = preCalcComponents(b,pi/4,pi/4,0,0);
%! X = createSineComponents(b,f,cM,0);
%! X2 = createSineComponents(b,f2,cM,0);
%! assert (rank(X2'*X2) == columns(X2))
%! assert (rank(X'*X) == columns(X))

