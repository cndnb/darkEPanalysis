function rtn = riceppf(cVal,sDev,pCN,cTol)
  xGuess = cVal + sDev;
  rtn = newtonMethod(xGuess,cVal,sDev,pCN,cTol);
endfunction

function rtn = newtonMethod(x,cVal,sDev,pCN,cTol)
  if(abs(1-marcumq(cVal/sDev,x/sDev) - pCN) < cTol)
  else
  %Newton's method, x_n+1 = x_n - f(x_n)/f'(x_n)
  newx = x - ((1-marcumq(cVal/sDev,x/sDev))/(besselj();
  endif
  rtn = x
endfunction
