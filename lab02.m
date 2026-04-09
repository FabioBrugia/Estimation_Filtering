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

