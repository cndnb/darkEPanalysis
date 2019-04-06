function rtn = riceppf(cVal,sDev,pCN,cTol)
	if(sDev <= 0)
		error("riceppf - sDev must be positive");
	endif
	if(cVal < 0)
		error("riceppf - cVal must be non-negative");
	endif
	xGuess = cVal + 2*sDev
	x = xGuess;

	count = 1;
	while(rCDF(x,cVal,sDev)/pCN > cTol & count < 100)
		if(x < 0)
			x=1/(10*sDev);
		endif
		x = x - ((rCDF(x,cVal,sDev) - pCN)/riceDistribution(x,cVal,sDev))
		fflush(stdout);
		count = count + 1;
	endwhile

	assert(1-marcumq(cVal/sDev,x/sDev),pCN,pCN*cTol);
	rtn = x;
endfunction

%Gives the rice probability density at each x for given central value and standard deviation
function rtn = riceDistribution(x,cVal,sDev)
	rtn = (x/(sDev)^2)*exp(-(x^2 + cVal^2)/(sDev^2))*besseli(0,(x*cVal)/(sDev^2));
endfunction

%Gives the cumulative probability of the rice distribution for x give cVal and sDev
function rtn = rCDF(x,cVal,sDev);
	rtn = 1 - marcumq(cVal/sDev,x/sDev);
endfunction

%!test
%! pCN = .95;
%! cTol = .01;
%! cVal = poissrnd(20); sDev = poissrnd(20) + 1;
%! xOut = riceppf(cVal,sDev,pCN,cTol);
%! assert(1-marcumq(cVal/sDev,xOut/sDev),pCN,cTol*pCN);
