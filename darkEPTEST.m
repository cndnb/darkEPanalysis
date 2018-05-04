t= 1:1e5; t=t';


tauT = [t,randn(length(t),1).*(1e-14)]; %Time domain torque (random noise)

%Pendulum constants
I = 378/1e7;                                                                    
f0 = 1.9338e-3;                                                                 
kappa = ((2*pi*f0)^2)*I;                                                        
Q = 1;                                                                     
T = 0; %No temp noise


preTheta = torqueSim(t,I,kappa,Q,T,tauT(:,2));

thetaT = [preTheta(:,1),preTheta(:,2)]; %Pulls out displacement values from torqueSim

tauW = psd(tauT(:,1),tauT(:,2)); %Frequency power spectrum of tau
thetaW = psd(thetaT(:,1),thetaT(:,2)); %Frequency power specturm of theta

figure(1);
plot(thetaT(:,1),thetaT(:,2),tauT(:,1),tauT(:,2),thetaT(:,1),preTheta(:,3)); %Time domain comparison
legend('Theta(t)','tau(t)','w');

figure(2);
loglog(tauW(:,1),tauW(:,2),thetaW(:,1),thetaW(:,2)); %Frequency domain comparison
 
figure(3);
tfVal = abs(transferFunction(tauW(:,1),kappa,f0,Q).^2); %ratio of power spectra vs power of transfer function
loglog(tauW(:,1),thetaW(:,2)./tauW(:,2),tauW(:,1),tfVal);

figure(4)
semilogx(tauW(:,1),thetaW(:,2)./(tauW(:,2).*tfVal));