function rtn = riceDistribution(xIn,cVal,sDev)
	if(nargin != 3)
		error("rtn = riceDistribution(xIn,cVal,sDev)");
	endif
	if(xIn < 0)
		error("riceDistribution - input value must be non-negative");
	endif
	if(cVal < 0)
		error("riceDistribution - central value must be non-negative");
	endif
	if(sDev < 0)
		error("riceDistribution - standard deviation must be positive");
	endif
	%Gives probability density for x given central value and standard deviation
	rtn = (xIn/(sDev)^2)*exp(-(xIn^2+cVal^2)/(sDev^2))*besseli(0,(xIn*cVal)/(sDev^2));
endfunction

%!test
%! t = (0:0.1:1000)';
%! area = ones(rows(t)-1,1);
%! for count = 1:rows(t)-1
%!  area(count,1) = .5*(riceDistribution(t(count),10,5) + riceDistribution(t(count + 1),10,5))*0.1;
%! endfor
%! assert(sum(area),1);
