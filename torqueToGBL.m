function rtn = torqueToGBL(torque)
	tFactor = (2e-11)*2*sqrt(3)./(0.04*0.02*9.8*0.037);
	rtn = torque.*tFactor;
endfunction
