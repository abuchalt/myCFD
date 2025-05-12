% ===================================================================
% The governing equation is DtQ + DxF = x in 0 < x < 1
% where Q = [U1;U2], Q = 0 at t = 0 and U1 = 0 at x = 0
%
% We perform a cell-centered FV discretization so that
%
% F1   Q2   F2   Q3   F3   Q4   F4   Q5   F5.....FN  QN+1  FN+1
% x----o----x----o----x----o----x----o----x......x----o----x
% Q1                                                       QN+2
%
% With this discretization, the cell-centers are defined as
% xc_i = (i - 1.5)*Delta_x, i = 2:N+1
%
% and similarly the cell-faces are defined as
% xf_i = (i - 1.0)*Delta_x, i = 1:N+1
%
% where Delta_x = 1.0/(N)
%
% According to the stencil used, the discretized set of ODEs become
%
% d(Q_i)/dt = S_i - (F_i - F_i-1)/Delta_x; i = 2:N+1
%
% where the fluxes are given as
%
% F_i = F((Q_i+1 + Q_i)/2.0); i = 2:N
% 
% and 
% F_1 = U_1 = 0 (BC at x = 0)
% d(QN+2)/dt = 1 - 2.0*(FN+1 - FN+1/2)/Delta_x (for Right Boundary)
%
% ===================================================================
%
clear all; clc;

% System Constants
gamma = 1.4;
M_in = 1.25;
M_out = 0.45;

% Boundary Conditions
P_in = 1;
rho_in = 1;
c_in = sqrt(gamma*P_in/rho_in);
u_in = M_in*c_in;

% Define variables
fprintf('Total number of cells = ')
N = input('');
x_f = zeros(1,N+1);
x_c = zeros(1,N+1);
Q   = zeros(3,N+2);
F   = zeros(3,N+1);
R   = zeros(3,N+2);
M = linspace(M_in, M_out, N+2);
rho_p1 = linspace(rho_in, 2*rho_in, N+2);
u_p1 = linspace(M_in*sqrt(gamma),M_out*sqrt(gamma), N+2);
P_p1 = rho_p1;
c_p1 = sqrt(gamma*P_p1./rho_p1);
Q(1,:) = rho_p1;
Q(2,:) = rho_p1.*u_p1;
Q(3,:) = P_p1/(gamma-1) + 0.5*rho_p1.*u_p1.^2;
dx = 10.0/N;
x_f(1,1:N+1) = 10.0*(0:N)/N;
x_c(1,2:N+1) = 10.0*((1:N)-0.5)/N;
x_c(1,N+2)   = 10.0;
%n_rk = 4;
%alfa = [1/4 1/3 1/2 1];
%cfl  = sqrt(8.0);
n_rk = 3;
alfa = [1/2 1/2 1];
cfl  = 1.95;
% n_rk = 2;
% alfa = [1/2 1];
% cfl  = 1.0;
epsilon = -15;

% Define Nozzle Area
for i = 1:N+2
    A(1,i) = 1.398 + 0.347*tanh(0.8*x_c(1,i)-3.2);
    r(1,i) = sqrt(A(1,i)/pi);
    dAdx(1,i) = 0.8*0.347*sech(0.8*x_c(1,i)-3.2)^2;
end

% File Info
mydir='C:\\Users\\Bucky\\Downloads\\FV2_Results';
subfolder='Ncell'+string(N)+'_CFL'+string(cfl);
mkdir(fullfile(mydir,subfolder));

dt   = 0.95*dx*cfl;
%
% Iterate until convergence
tTot = 0;
niter = 1;
Residual(niter) = 0;
iteration(niter) = 1;
% plot(x_c, rho_p1);
% draw now;
% pause;
while Residual(niter) > epsilon
    tStart = tic;
    niter = niter + 1;
    iteration(niter) = niter;
% for niter = 1:300
    %
    % RK integration
    Q0 = Q;
    for k = 1:n_rk
        dt = 10;
        for i = 1:N+2
            rho_p1(1,i) = Q(1,i);
            u_p1(1,i) = Q(2,i)/Q(1,i);
            P_p1(1,i) = (gamma-1)*(Q(3,i)-0.5*(Q(2,i)^2/Q(1,i)));
            c_p1(1,i) = sqrt(gamma*P_p1(1,i)/rho_p1(1,i));
            dt = min(dt,0.95*dx*cfl/(abs(u_p1(1,i))+c_p1(1,i)));
        end
        %
        % Compute Conservation Variables at Faces
        for i = 2:N
            Qface(1,i) = 0.5*(Q(1,i)+Q(1,i+1));
            Qface(2,i) = 0.5*(Q(2,i)+Q(2,i+1));
            Qface(3,i) = 0.5*(Q(3,i)+Q(3,i+1));
        end
        % And at Boundaries (just remember dx is half-size for these)
        Qface(1,1) = Q(1,1);
        Qface(2,1) = Q(2,1);
        Qface(3,1) = Q(3,1);
        Qface(1,N+1) = Q(1,N+1);
        Qface(2,N+1) = Q(2,N+1);
        Qface(3,N+1) = Q(3,N+1);
        %
        % Compute Fluxes
        % for i = 2:N
        for i = 1:N+1
            F(1,i) = Qface(2,i);
            F(2,i) = (Qface(2,i)^2/Qface(1,i))*(3-gamma)/2 + (gamma-1)*Qface(3,i);
            F(3,i) = (gamma*Qface(2,i)*Qface(3,i)/Qface(1,i)) - ((gamma-1)/2)*Qface(2,i)^3/Qface(1,i)^2;
        end
        %
        % Compute Source Terms
        S(:,1:N+2) = 0.0;
        for i = 1:N+2 
            S(1,i) = -Q(2,i)*(dAdx(1,i)/A(1,i));
            S(2,i) = -(Q(2,i)^2/Q(1,i))*(dAdx(1,i)/A(1,i));
            S(3,i) = -(Q(3,i)/Q(1,i) + (gamma-1)*(Q(3,i)-0.5*Q(2,i)^2/Q(1,i)))*(dAdx(1,i)/A(1,i));
        end
        %
        % Compute Residuals
        R = 0.0;
        for i = 2:N+1
            R(1,i) = S(1,i) - (F(1,i) - F(1,i-1))/dx;
            R(2,i) = S(2,i) - (F(2,i) - F(2,i-1))/dx;
            R(3,i) = S(3,i) - (F(3,i) - F(3,i-1))/dx;
        end
        %
        % Take care of the Residuals at the ends for Q1 and Q_N+2
        F_right(1,1) = Q(2,2);
        F_right(2,1) = (Q(2,2)^2/Q(1,2))*(3-gamma)/2 + (gamma-1)*Q(3,2);
        F_right(3,1) = (gamma*Q(2,2)*Q(3,2)/Q(1,2)) - ((gamma-1)/2)*Q(2,2)^3/Q(1,2)^2;
        R(:,1) = S(:,1) - (F_right(:,1)-F(:,1))/(dx/2);

        F_left(1,1) = Q(2,N+1);
        F_left(2,1) = (Q(2,N+1)^2/Q(1,N+1))*(3-gamma)/2 + (gamma-1)*Q(3,N+1);
        F_left(3,1) = (gamma*Q(2,N+1)*Q(3,N+1)/Q(1,N+1)) - ((gamma-1)/2)*Q(2,N+1)^3/Q(1,N+1)^2;
        R(:,N+2) = S(:,N+2) - (F(:,N+1)-F_left(:,1))/(dx/2);
        %
        % Update Solution
        Q = Q0 + dt*alfa(k)*R;
        %
        % Apply BC - Use these for now! I'll go over the derivation
        % 
        % Supersonic Inlet
        % x = 0.0 BC
        Q(1,1) = rho_in;
        Q(2,1) = rho_in*u_in;
        Q(3,1) = P_in/(gamma-1) + 0.5*rho_in*u_in^2;
        %
        % Subsonic Outlet
        % x = 10.0 BC
        % Predicted Values
        rho_p = Q(1,N+2);
        u_p = Q(2,N+2)/Q(1,N+2);
        P_p = (gamma-1)*(Q(3,N+2)-0.5*(Q(2,N+2)^2/Q(1,N+2)));
        c_p = sqrt(gamma*P_p/rho_p);
        % corrected
        c = (u_p+(2/(gamma-1))*c_p)/(M_out+(2/(gamma-1)));
        rho = rho_p*(c/c_p)^(2/(gamma-1));
        P = c^2*rho/gamma;
        u = M_out*c;
        % sub back into conservation vars
        Q(1,N+2) = rho;
        Q(2,N+2) = rho*u;
        Q(3,N+2) = P/(gamma-1)+0.5*rho*u^2;
    end
    tTot = tTot + toc(tStart);
    %
    Residual(niter) = log10(norm(R)/(2*(N+2)));           
    % figure(1);
    % plot(x_c,M);
    % hold on;
    % plot(x_c,r);
    % legend('velocity', 'nozzle profile')
    % axis([0 1 0 0.6]);
    % ylabel('u');
    % xlabel('x');
    % title('Nozzle');
    figure(2);
    plot(x_c,Q(1,:),x_c,Q(2,:),x_c,Q(3,:))
    axis([0 10 0 4]);
    drawnow;
    pause(0.01);
    fprintf(1,'Iteration Number = %g Residual = %g\n',niter,Residual(niter));
end

figure(2);
plot(iteration, Residual, 'bo');
grid on;
ylabel('log10(residual)');
xlabel('iteration');
title('Convergence Behavior');

saveas(figure(1),fullfile(mydir,subfolder,subfolder+'_soln.jpg'));
saveas(figure(2),fullfile(mydir,subfolder,subfolder+'_conv.jpg'));

% And Iteration Info
fid = fopen(fullfile(mydir,subfolder,'iter.txt'),'wt');
fprintf(fid, 'Total Iterations: %s \nTotal CPU-time: %s s', string(niter), string(tTot));
fclose(fid);