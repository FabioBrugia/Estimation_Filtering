 clc
clear all
close all

s = open("sensor.mat");
z = s.z;
Vz = s.Vz;

%%
figure, plot(z,Vz,'o'),
axis([min(z),max(z),-10,10]), grid, title('Position transducer'),
xlabel('Position{\it z} in m'), ylabel('Voltage{\it V_z} in V')
%%
i1 = 11; 
i2 = 35;
z_lin = z(i1:i2);
Vz_lin = Vz(i1:i2); 
N_lin = length(z_lin);
Phi = [Vz_lin, ones(i2-i1+1,1)];

p = Phi\z_lin;

Kt=1/p(1)
Vo=-p(2)/p(1)

Vz0=linspace(-10,10,1000);
z_hat=p(1)*Vz0+p(2);
figure, plot(z,Vz,'o', z_hat,Vz0,'-', ...
z_lin(1)*[1,1],[-10,10],'-', z_lin(end)*[1,1],[-10,10],'-'),
axis([min(z),max(z),-10,10]), grid, title('Position transducer'),
ylabel('Voltage{\it V_z} in V'),
text(0.013,-9,'\fontsize{20} linearity interval')
%%
%% Computation of the parameter confidence intervals (noise variance derived from priors)
% Step 5: computation of the confidence intervals
% Definition of the x% confidence intervals:
% x=95.4 => k=2 ("2 sigma"); x=99.7 => k=3 ("3 sigma")
k_e=2; k_p=2;
noise_max=5e-4;
sigma_e=noise_max/k_e
Sigma_e=sigma_e^2*eye(N_lin);
Sigma_p=inv(Phi'*inv(Sigma_e)*Phi);
sigma_p=sqrt(diag(Sigma_p));
delta_p=k_p*sigma_p;
p_min=p-delta_p;
p_max=p+delta_p;
Kt_min=1/p_max(1)
Kt_max=1/p_min(1)
Vo_min=-p_max(2)/p_min(1)
Vo_max=-p_min(2)/p_max(1)
% Step 6: graphical comparison of the results
z_min=min([p_max(1)*Vz0; p_min(1)*Vz0])+p_min(2);
z_max=max([p_max(1)*Vz0; p_min(1)*Vz0])+p_max(2);
figure, plot(z,Vz,'o', z_hat,Vz0,'-', ...
z_min,Vz0,'--', z_max,Vz0,'-.', ...
z_lin(1)*[1,1],[-10,10],'-',z_lin(end)*[1,1],[-10,10],'-')
axis([min(z),max(z),-10,10]), grid,
title('Position transducer (noise variance derived from priors)'),
xlabel('Position{\it z} in m'), ylabel('Voltage{\it V_z} in V')
%% Computation of the parameter confidence intervals (noise variance estimated from data)
% Step 7: computation of the confidence intervals
sigma_e_ml=sqrt((z_lin-Phi*p)'*(z_lin-Phi*p)/N_lin)
Sigma_e=sigma_e_ml^2*eye(N_lin);
Sigma_p=inv(Phi'*inv(Sigma_e)*Phi);
sigma_p=sqrt(diag(Sigma_p));
delta_p=k_p*sigma_p;
p_min=p-delta_p;
p_max=p+delta_p;
Kt_min=1/p_max(1)
Kt_max=1/p_min(1)
Vo_min=-p_max(2)/p_min(1)
Vo_max=-p_min(2)/p_max(1)
% Step 8: graphical comparison of the results
z_min=min([p_max(1)*Vz0; p_min(1)*Vz0])+p_min(2);
z_max=max([p_max(1)*Vz0; p_min(1)*Vz0])+p_max(2);
figure, plot(z,Vz,'o', z_hat,Vz0,'-', ...
z_min,Vz0,'--', z_max,Vz0,'-.', ...
z_lin(1)*[1,1],[-10,10],'-',z_lin(end)*[1,1],[-10,10],'-')
axis([min(z),max(z),-10,10]), grid,
8
title('Position transducer (noise variance estimated from data)'),
xlabel('Position{\it z} in m'), ylabel('Voltage{\it V_z} in V'),


%%
% a
L = [Vz.^3 Vz.^2 Vz Vz.^0];
p3a = L\z;
z_pol = polyval(p3a,Vz0);
% b
L=[z.^3, z.^2, z, z.^0];
p3b=L\Vz
z0=linspace(min(z),max(z),1000);
Vz_pol=polyval(p3b,z0);

figure, plot(z,Vz,'o',z_pol,Vz0,'--',z0,Vz_pol,'-'),
axis([min(z),max(z),-10,10]), grid, title('Position transducer'),
xlabel('Position{\it z} in m'), ylabel('Voltage{\it V_z} in V'),
print -deps sensore4.eps


%% -------------------------------------------------------
% 3rd order polymonial model by means of polyfit
%-------------------------------------------------------
modello_ordine3=polyfit(z,Vz,3);
Vz_modello_ordine3=polyval(modello_ordine3,z0);

