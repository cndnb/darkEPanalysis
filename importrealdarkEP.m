clear d;

%imports data
d = load('run0160.dat');
t = (1:rows(d))';
%d(:,1) = 86400.*d(:,1); %Convert from Julian days to seconds
%d(:,1) = d(:,1) - d(1,1).*ones(rows(d),1);
d(:,4) = (14e-6).*d(:,4); %Convert from pixels to radians
d = [t,d(:,4),d(:,8)];
