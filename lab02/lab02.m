clc
clear all
close all

load sensor_lab02.mat
figure, plot(z,Vz,'o'),
axis([min(z),max(z),-10,10]), grid, title('Position transducer'),
xlabel('Position{\it z} in m'), ylabel('Voltage{\it V_z} in V')

%%

i1 = 14;
i2 = 36;

z_lin = z(i1:i2);
Vz_lin = Vz(i1:i2);

N_lin = length(z_lin);

Phi = [Vz_lin, ones(N_lin,1)];

p = Phi\z_lin

A=pinv(Phi); % more realiable than inv(Phi'*Phi)*Phi'
p_=A*z_lin; % Form #2: using the pseudoinverse matrix
Kt_=1/p_(1)
Vo_=-p_(2)/p_(1)
% Step 4: graphical comparison of the results
Vz0=linspace(-10,10,1000);
z_hat=p(1)*Vz0+p(2);
figure, plot(z,Vz,'o', z_hat,Vz0,'-', ...
z_lin(1)*[1,1],[-10,10],'-', z_lin(end)*[1,1],[-10,10],'-'),
axis([min(z),max(z),-10,10]), grid, title('Position transducer'),
text(0.015,-9,'\fontsize{20} linearity interval')
%%
% bound on the greatest measurement error of the position z
format compact, format short e
epsilon=0.5e-3
% EUI algorithm in compact form
for ind=1:length(p),
    p_min(ind,1)=A(ind,:)*(z_lin-epsilon*sign(A(ind,:))');
end
p_max=2*p-p_min;
EUI=[p_min, p_max]
% EUI algorithm in extended form, using a loop (equivalent)
p_min_=zeros(size(p)); ; p_max_=zeros(size(p));
for k=1:length(A),
    p_min_(1)=p_min_(1)+A(1,k)*(z_lin(k)-epsilon*sign(A(1,k)));
    p_min_(2)=p_min_(2)+A(2,k)*(z_lin(k)-epsilon*sign(A(2,k)));
    p_max_(1)=p_max_(1)+A(1,k)*(z_lin(k)+epsilon*sign(A(1,k)));
    p_max_(2)=p_max_(2)+A(2,k)*(z_lin(k)+epsilon*sign(A(2,k)));
end
EUI_=[p_min_, p_max_]
format compact, format short
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
title('Position transducer (uncertainty region due to EUI)'),
xlabel('Position{\it z} in m'), ylabel('Voltage{\it V_z} in V'),

%% Computation of the PUIs
options_old=optimset('linprog'); % default value for 'Algorithm' option is 'interior-point'
options_new=optimset(options_old, 'Algorithm','interior-point');

PUI(1,1)=[1; 0]'*linprog( [1; 0],[Phi; -Phi],[z_lin+epsilon; -z_lin+epsilon],[],[],[],[],options_new);
PUI(2,1)=[0; 1]'*linprog( [0; 1],[Phi; -Phi],[z_lin+epsilon; -z_lin+epsilon],[],[],[],[],options_new);
PUI(1,2)=[1; 0]'*linprog(-[1; 0],[Phi; -Phi],[z_lin+epsilon; -z_lin+epsilon],[],[],[],[],options_new);
PUI(2,2)=[0; 1]'*linprog(-[0; 1],[Phi; -Phi],[z_lin+epsilon; -z_lin+epsilon],[],[],[],[],options_new);


format compact, format short e
PUI
format compact, format short
Kt_min_PUI=1/PUI(1,2)
Kt_max_PUI=1/PUI(1,1)
Vo_min_PUI=-PUI(2,2)/PUI(1,1)
Vo_max_PUI=-PUI(2,1)/PUI(1,2)

%%step 8

format compact, format short e
p_central = [(PUI(1,1) + PUI(1,2))/2; (PUI(2,1)+PUI(2))/2]

format compact, format short
kt_central = 1/p_central(1)
Vo_central = -p_central(2)/p_central(1)

%% step 9 comparisione of the result 

z_hat_central=p_central(1)*Vz0+p_central(2);
z_min_central=min([PUI(1,2)*Vz0; PUI(1,1)*Vz0])+PUI(2,1);
z_max_central=max([PUI(1,2)*Vz0; PUI(1,1)*Vz0])+PUI(2,2);
figure, plot(z,Vz,'o', z_hat_central,Vz0,'-', ...
z_min_central,Vz0,'--', z_max_central,Vz0,'-.', ...
z_lin(1)*[1,1],[-10,10],'-',z_lin(end)*[1,1],[-10,10],'-')
axis([min(z),max(z),-10,10]), grid,
title('Position transducer (uncertainty region due to PUI)'),
xlabel('Position{\it z} in m'), ylabel('Voltage{\it V_z} in V'),
legend('Experimental data (z, V_z)_{ }^{ }', 'z \approx \theta_{1}^{C} * V_{z} + \theta_2^{C}','z \approx min(\theta_{1}^{M} * V_{z} , \theta_{1}^{m} * V_{z}) + \theta_{2}^{m}','z \approx max(\theta_{1}^{M} * V_{z} , \theta_{1}^{m} * V_{z}) + \theta_{2}^{M}')
text(0.015,-9,'\fontsize{20} linearity interval')
