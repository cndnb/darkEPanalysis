function ret = removeHighError(t,y,err,threshold)
  removeCheck = zeros(length(t),1);
  count = 0;
  for n = 1:length(t)
    if (err(n) > threshold)
      removeCheck(n,1) = 1;
      count = count + 1;
    endif
  endfor
  results = ones(length(t)-count,3);
  count = 0;
  for i = 1:length(t)
    if (removeCheck(i,1)==1)
      count = count+1;
    else
       results(i,:) = [t(i+count),y(i+count),err(i+count)];
    endif
  endfor
  
  ret = results;
endfunction