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
