function derivative = symmetricDerivative(data)
  dd = diff(data);

	derivative = ( dd(1:end-1,:) + dd(2:end,:) ) / 2;
  
endfunction

%!test
%!
%! t = 8*(transpose(1:100));
%! d = symmetricDiff(t);
%! assert( d == 8*ones(98,1))

%!test
%!
%! %Checking against sinusoids
%! t = transpose(1:10000);
%! f = 0.0011;
%! omega = 2 * pi * f;
%! s = sin( omega * t );
%! c = omega * cos( omega * t );
%! d = symmetricDiff(s);
%! plot( (d-c(2:end-1)) ./ d )
%! assert( (d - c(2:end-1) ) ./d, zeros(size(d)), 1e-5 )

%Old code
%out = ones(length(data),1);
%out(1,1) = data(2,2)-data(2,1);
%out(length(data),1) = data(length(data),2)-data(length(data)-1,2);
%for i = 2:length(data)-1
%  out(i,1) = (data(i+1,2)-data(i-1,2))/2;
%endfor
%derivative = [data(:,1), out];