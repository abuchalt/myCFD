clear all; clc;
%
% System Constants
gamma = 1.4;
M_in = 1.25;
M_out = 0.45;
N=500;
%
% File Info
mydir='C:\\Users\\Bucky\\Downloads\\FV2_Results';
subfolder='Ncell'+string(N)+'Jameson';
%
% Domain
x_c = zeros(1,N+1);
dx = 10.0/N; % Domain Defined over 10.0 untis
x_c(1,2:N+1) = 10.0*((1:N)-0.5)/N;
x_c(1,N+2)   = 10.0;
%
% Define Nozzle Area
for i = 1:N+2
    A(1,i) = 1.398 + 0.347*tanh(0.8*x_c(1,i)-3.2);
    r(1,i) = sqrt(A(1,i)/pi);
    dAdx(1,i) = 0.8*0.347*sech(0.8*x_c(1,i)-3.2)^2;
end
%
% Solution
myQ = matfile(fullfile(mydir,subfolder,subfolder+'_Q.mat'));
Q = myQ.Q;
QSmooth = Q; % Conservation State Variables
%
% Second-Order Smoothing Filter
for i = 2:N+1
    QSmooth(:,i) = 0.25*(Q(:,i-1)+2*Q(:,i)+Q(:,i+1));
end
Q = QSmooth;
%
% Re-Plot  
figure(1);
plot(x_c,Q(1,:))
hold on;
plot(x_c,Q(2,:))
hold on;
plot(x_c,Q(3,:))
legend('Q1', 'Q2', 'Q3')
axis([0 10 0 6]);
xlabel('x');
title('Smoothed Conservation Variables');
drawnow;
hold off;
%
saveas(figure(1),fullfile(mydir,subfolder,subfolder+'_smoothedCons.jpg'));
%
for i = 1:N+2
    rho_p1(1,i) = Q(1,i);
    u_p1(1,i) = Q(2,i)/Q(1,i);
    P_p1(1,i) = (gamma-1)*(Q(3,i)-0.5*(Q(2,i)^2/Q(1,i)));
    c_p1(1,i) = sqrt(gamma*P_p1(1,i)/rho_p1(1,i));
end
M = u_p1./c_p1;
for i = 1:N+2
    c(1,i) = (u_p1(1,i)+(2/(gamma-1))*c_p1(1,i))/(M(1,i)+(2/(gamma-1)));
    rho(1,i) = rho_p1(1,i)*(c(1,i)/c_p1(1,i))^(2/(gamma-1));
    P(1,i) = c(1,i)^2*rho(1,i)/gamma;
    u(1,i) = M(1,i)*c(1,i);
end
%
figure(2);
plot(x_c,M);
hold on;
plot(x_c,r);
legend('mach number', 'nozzle profile')
axis([0 10 0 2.5]);
xlabel('x');
title('Smoothed Flow Through a Nozzle');
drawnow;
hold off;
%
saveas(figure(2),fullfile(mydir,subfolder,subfolder+'_smoothedSoln2.jpg'));