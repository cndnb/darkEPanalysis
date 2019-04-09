function rtn = riceppf(cVal,sDev,pCN,cTol,xGuess)
	if(sDev <= 0)
		error("riceppf - sDev must be positive");
	endif
	if(cVal < 0)
		error("riceppf - cVal must be non-negative");
	endif

	%Prevents running forever
	maxCount = 100;
	
	%Starts method with guess value
	x = xGuess;

	%initialize variables
	count = 1;
	newX = 0;
	while(count < maxCount)
		%Sets new value with Newton's method
		newX = x - ((rCDF(x,cVal,sDev) - pCN)/riceDistribution(x,cVal,sDev));
		%Prevents overshooting from breaking the program
		if(newX < 0)
			x=1/(10*sDev);
		else
			%Recursion step
			x = newX;
		endif
		%If the value gives accurate enough result, we are done
		if(abs((rCDF(x,cVal,sDev)/pCN) - 1) < cTol)
			break;
		endif
		count = count + 1;
	endwhile
	
	%If count == maxCount the return value does not necessarily have percent accuracy of cTol
	if(count == maxCount)
		warning("riceppf - precision not accurate");
	endif
	
	%Returns value
	rtn = x;
endfunction

%!test
%! pCN = .95;
%! cTol = eps
%! cVal = poissrnd(40) 
%! sDev = poissrnd(10) + 1
%! xOut = riceppf(cVal,sDev,pCN,cTol,cVal + 2*sDev)
%! assert(rCDF(xOut,cVal,sDev),pCN,cTol*pCN);
