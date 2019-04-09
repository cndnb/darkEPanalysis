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
	rtn = (xIn./(sDev).^2).*exp(-(xIn.^2+cVal.^2)./(2*sDev.^2)).*besseli(0,(xIn.*cVal)/(sDev.^2));
endfunction

%!test
%! interval = 0.01;
%! cVal = poissrnd(20);
%! sDev = poissrnd(10) + 1;
%! t = (0:interval:cVal + 10*sDev)';
%! area = ones(rows(t)-1,1);
%! for count = 1:rows(t)-1
%!  area(count,1) = riceDistribution((t(count) + t(count+1))./2,cVal,sDev)*interval;
%! endfor
%! assert(sum(area),1,1e-5);
