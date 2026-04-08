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

phi = [Vz_lin, ones(i2-i1+1,1)];

theta = phi\z_lin;

Kt=1/theta(1)
Vo=-theta(2)/theta(1)

Vz0=linspace(-10,10,1000);
z_hat=theta(1)*Vz0+theta(2);
figure, plot(z,Vz,'o', z_hat,Vz0,'-', ...
z_lin(1)*[1,1],[-10,10],'-', z_lin(end)*[1,1],[-10,10],'-'),
axis([min(z),max(z),-10,10]), grid, title('Position transducer'),
ylabel('Voltage{\it V_z} in V'),
text(0.013,-9,'\fontsize{20} linearity interval')
%%





