function rtn = preCalcComponents(timeData,seattleLat,seattleLong,compassDir,startTime)
	if (nargin != 5)
		error("rtn = preCalcComponents(time,Lat,Long,compassDir,startTime)");
	endif
	
	%Creates array to store precalculated multiples for final component calculation
	%Z = 1, perpX = 2, paraX = 3
	constants = ones(rows(timeData),3);
	
        omegaEarth = 2*pi*(1/86164.0916);
	%Precalculates large argument and takes mod so that sin/cos calculate it more accurately
	latAngle = omegaEarth.*timeData + seattleLong + omegaEarth*startTime;
	fSLAT = seattleLat;
	fSLONG = seattleLong;
	%compass direction is already in the range 0 to 2pi

	%Z constant
	constants(:,1) = cos(fSLAT)*cos(compassDir).*ones(rows(timeData),1);

	%PerpX constant
	constants(:,2) = (-sin(compassDir).*cos(latAngle)) .- (sin(fSLAT)*cos(compassDir).*sin(latAngle));

	%ParaX constant
	constants(:,3) = (sin(compassDir).*sin(latAngle)) .- (sin(fSLAT)*cos(compassDir).*cos(latAngle));

	%Returns constant array
	rtn = constants;
endfunction
	
%!test
%! omegaEarth = 2*pi*(1/86164.0916);
%! cM = preCalcComponents([0;pi/(2*omegaEarth)],0,0,pi/4,0);
%! assert(cM,[sqrt(2)/2,-sqrt(2)/2,0;sqrt(2)/2,0,sqrt(2)/2],eps);	

%!test
%! assert(preCalcComponents(0,0,0,0,0) == [1,0,0]);

%!test
%! assert(preCalcComponents(0,0,0,-pi/2,0),[0,1,0],eps);

%!test
%! assert(preCalcComponents(0,0,pi/2,pi/2,0),[0,0,1],eps);
