function ret = dayDivision(data,daysInclude,dayLength)
  if(nargin != 3)
    error('cellArray = dayDivision(data,daysInclude,dayLength)');
  endif
  
  %Finds index of day rollover
  modCount = mod(data(:,1) - data(1,1).*ones(rows(data),1),dayLength);
  modCount = [modCount(1,1);[modCount(2:end,1) - modCount(1:end-1,1) + ones(rows(modCount)-1,1)]];
  indVal = find(modCount./abs(modCount) - 1);
  %Sets initial value so that first day always starts at the beginning of data
  if(indVal(end,1) != rows(data))
    indVal = [indVal;rows(data)];
  endif
  %Remove starting index from maxDays
  maxDays = rows(indVal) - 1;
  
  indVal
  fflush(stdout);
  
  %Initialize accumulation array
  dataDivisions = cell(maxDays,1);
  
  %Assigns each interval to a cell
  for count = 1:rows(indVal) - 1
    dataDivisions{count,1} = data(indVal(count):indVal(count + 1) - 1,:);
  endfor
  
  %Return
  numDays = 0;
  if (daysInclude == 0) %Fits the maximum number of whole days in the data
    numDays = maxDays;
  elseif (daysInclude > maxDays)
    numDays = maxDays;
  else
    numDays = daysInclude;
  endif
  ret = dataDivisions(1:numDays,1);

endfunction

%!test
%! t = 1:10000; t=t'; %Create time
%! fData = [t,randn(rows(t),1)]; %Create fake data [t,Y]
%! dayLength = 100;
%! divData = dayDivision(fData,0,dayLength); %0 => Include everything
%! checkData = cell2mat(divData(:,1));
%! assert(checkData,fData(1:rows(checkData),:));
