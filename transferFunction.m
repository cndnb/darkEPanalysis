function ret = transferFunction(inFreq,kappa,f0,Q)
  if (columns(inFreq)!= 1)
    error("transferFunction can only take a frequency column vector");
  endif
  if (nargin != 4)
    usage("tau(frequency) = transferFunction(frequency,kappa,resonanceFreq,Q)");
  endif
  
  ret = (1/kappa)./((1 .- (inFreq./f0).^2).+((sqrt(-1).*inFreq)./(Q*f0)));
endfunction

%!test
%! inFreq = ones(10,1);
%! assert(abs(transferFunction(inFreq,1,2,2) - ((6-2*sqrt(-1))/5).*ones(10,1))<eps)

%!test
%! inFreq = [1/10,1/100,1/100]'; %Take some frequencies
%! kappa = pi; f0 = e; Q = sqrt(2);
%! fVal = transferFunction(inFreq,kappa,f0,Q);
%! rVal = ones(rows(inFreq),1);
%! for ind = 1:rows(inFreq)
%!  rVal(ind) = (1/kappa)/((1-(inFreq(ind)/f0)^2)+((sqrt(-1)*inFreq(ind))/(Q*f0)));
%! endfor
%! for count = 1:rows(inFreq)
%!  assert(fVal(count) == rVal(count)) %Makes sure each column is the right number
%! endfor

%!test
%! t = 1:1e5;t=t';
%! tw = 5;
%! O = torqueSim(t,I,kappa,Q,T,tw)