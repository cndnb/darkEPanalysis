function ret = twoPointTransfer(inFreq,f0) %Update to work with discrete data
  if (nargin != 2)
    usage("tau(frequency) = twoPointTransfer(frequency,resonanceFreq)");
  endif
  if (columns(inFreq)!= 1)
    error("transferFunction can only take a frequency column vector");
  endif
  
  ret = cos((pi.*inFreq)/(2*f0));
endfunction

%!test
%! freq = ones(10,1);
%! f0 = randn;
%! eVal = cos(pi/(2*f0));
%! cVal = twoPointTransfer(freq,f0);
%! assert(cVal,eVal.*ones(rows(freq),1));

%!test
%! freqRange = [1/100,1/10,1]';
%! f0 = 1.9338e-3;
%! eVal = cos((pi.*freqRange)/(2*f0));
%! cVal = twoPointTransfer(freqRange,f0);
%! assert(cVal,eVal);
