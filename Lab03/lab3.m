clc
clear all
close all

load resistor_data_1.mat
N = length(i);
%%
figure,
plot(i,V)
axis([min(i),max(i),0,35]), grid, title('Measured data and approximations'),
xlabel('Current{\it i} in A'), ylabel('Voltage{\it V} in V'),

%%

R_SM = mean(V./i);

V_SM = R_SM*i;

%%
figure,
plot(i,V,'r',i,V_SM,'-b')
axis([min(i),max(i),0,35]), grid, title('Measured data and approximations'),
xlabel('Current{\it i} in A'), ylabel('Voltage{\it V} in V'),

%% Estimation Using LS
Phi = i;
R_LS = Phi\V
A = pinv(Phi);
sigma_e = 5;
sigma_ee = sigma_e^2*eye(N);
sigma_R_LS = pinv(Phi)*sigma_ee*pinv(Phi)';
%%
V_LS = R_LS*i;
figure, plot(i,V,'og', i,V_SM,'-b', i,V_LS,'-r'),
axis([min(i),max(i),0,35]), grid, title('Measured data and approximations'),
xlabel('Current{\it i} in A'), ylabel('Voltage{\it V} in V'),
%% Gaussian

R_GM = inv(Phi'*inv(sigma_ee)*Phi)*Phi'*inv(sigma_ee)*V
Sigma_R_GM = inv(Phi'*inv(sigma_ee)*Phi)

V_GM = R_GM*i;
figure, plot(i,V,'ok', i,V_SM,'-b', i,V_LS,'-r', i,V_GM,'--k'),
axis([min(i),max(i),0,35]), grid, title('Measured data and approximations'),
xlabel('Current{\it i} in A'), ylabel('Voltage{\it V} in V')
%% Bysian method


R_bar = R_GM; % Try also with R_GM
Sigma_RR = [10, 1, 0.1, 0.01];
V_bar = R_bar*i;
for k=1:length(Sigma_RR),
    Sigma_VV = Sigma_RR(k)* i*i' + sigma_ee;
    Sigma_RV = Sigma_RR(k)*i';
    R_B(k) = R_bar + Sigma_RV * inv(Sigma_VV) * (V-V_bar);
    Sigma_R_B(k) = Sigma_RV * inv(Sigma_VV) * Sigma_RV';
    Sigma_R_minus_R_B(k) = Sigma_RR(k) - Sigma_R_B(k);
    V_B(:,k) = R_B(k)*i;
end
Sigma_RR, R_B, Sigma_R_B, Sigma_R_minus_R_B

%%
    