function rtn = rCDF(xIn,cVal,sDev)
	rtn = 1 - marcumq(cVal/sDev,xIn/sDev);
endfunction

%!test
%! assert(rCDF(1000,1,1), 1);

%!test
%! assert(rCDF(eps,1,1),0);
