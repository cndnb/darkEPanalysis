gZ = torqueToGBL(FINALAMP(:,2));
gPerpX = torqueToGBL(FINALAMP(:,3));
gParaX = torqueToGBL(FINALAMP(:,4));
errZ = torqueToGBL(FINALERR(:,2));
errPerpX = torqueToGBL(FINALERR(:,3));
errParaX = torqueToGBL(FINALERR(:,4));

f = 1:10^4; f=f./10^5;
f=f';
g1 = (2e-24).*ones(rows(f),1);
g2 = (1.5e-25).*ones(rows(f),1);
g3 = (9e-27).*ones(rows(f),1);

figure(1);
%loglog(FINALAMP(:,1),[gZ,gZ - errZ,gZ+errZ,gLim],f,[g1,g2,g3]);
loglog(FINALAMP(:,1),[gZ,gLim],f,[g1,g2,g3]);
legend('Torsion Filter Fit','Thermal Limit','static EP tests','reanalysis','next run');
title('Exclusion Plot for Z coupling');
%legend('Resonance Fit','err+','err-','Thermal Limit','static EP tests','reanalysis','next run');
xlabel('Frequency (Hz)');
ylabel('g_{B-L}');

figure(2);
loglog(FINALAMP(:,1),[gPerpX,gLim],f,[g1,g2,g3]);
title('Exclusion Plot for Perpendicular to Gamma coupling');
legend('Torsion Filtered','Thermal Limit','static EP tests','reanalysis','next run');
xlabel('Frequency (Hz)');
ylabel('g_{B-L}');

%figure(3);
%loglog(FINALAMP(:,1),[gParaX,gParaX2,gLim],f,[g1,g2,g3]);
%title('Exclusion Plot for Parallel to Gamma coupling');
%legend('Resonance Fit','Torsion Filtered','Thermal Limit','static EP tests','reanalysis','next run');
%xlabel('Frequency (Hz)');
%ylabel('g_{B-L}');
