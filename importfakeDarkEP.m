clear d;

%imports data
d = load('fakeDarkEPAugust92017.dat');
d = [d,ones(rows(d),1)];