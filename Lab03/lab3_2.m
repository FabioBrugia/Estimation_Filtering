clc, clear all, close all

load AR_data.mat

L = length(y);
%% 
figure,
plot(1:L, y,'-b')
axis([0,L,0,100]), grid, title('Measured data and approximation'),
xlabel('Time{\it t}'), ylabel('Output{\it y}'),
legend('Measured data')

%% step 3 computation of the estimated paramater

Y = y(3:L);
Phi = [-y(2:L-1), -y(1:L-2), ones(L-2,1)];

theta = Phi\Y;

a1=theta(1)
a2=theta(2)
b=theta(3)

%%
sigma_e=sqrt(40);
Sigma_theta=sigma_e^2*inv(Phi'*Phi) % Variance matrix of the estimates
sigma_theta=sqrt(diag(Sigma_theta)) % Standard deviations of the estimates
%%
y_LS = Phi*theta;
figure, plot(1:L,y,'-b', 3:L,y_LS,'-r'),
axis([0,L,0,100]), grid, title('Measured data and approximation'),
xlabel('Time{\it t}'), ylabel('Output{\it y}'),
legend('Measured data','Least Squares estimate')

