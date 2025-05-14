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
%
% System Constants
gamma = 1.4;
M_in = 1.25;
M_out = 0.45;
%
% Mac-Baldwin Dissipation Factor
epsilon = 0.8;
%
% Jameson Dissipation Factors
alpha2 = 1/4;
alpha4 = 1/256;
%
% Boundary Conditions
P_in = 1;
rho_in = gamma;
c_in = sqrt(gamma*P_in/rho_in);
u_in = M_in*c_in;
%
% Define variables
fprintf('Total number of cells = ')
N = input('');
%
x_f = zeros(1,N+1); % Face domain
x_c = zeros(1,N+1); % Center domain
Q   = zeros(3,N+2); % Conservation State Variables
F   = zeros(3,N+1); % Fluxes
R   = zeros(3,N+2); % Residuals
%
M = linspace(M_in, M_out, N+2); % Interpolate Mach IG
u_p1 = linspace(M_in, M_out, N+2); % Interpolate Velocity IG (c hard coded to 1)
P_p1 = linspace(1, 1, N+2); % Interpolate Pressure IG
c_p1 = M./u_p1; % Calculate c from interp IG
rho_p1 = gamma*P_p1./(c_p1.^2); % Calculate Density IG
%
Q(1,:) = rho_p1; % Build IG
Q(2,:) = rho_p1.*u_p1;
Q(3,:) = P_p1/(gamma-1) + 0.5*rho_p1.*u_p1.^2;
%
dx = 10.0/N; % Domain Defined over 10.0 untis
x_f(1,1:N+1) = 10.0*(0:N)/N;
x_c(1,2:N+1) = 10.0*((1:N)-0.5)/N;
x_c(1,N+2)   = 10.0;
%
% RK Params
n_rk = 4;
alfa = [1/4 1/3 1/2 1];
cfl  = sqrt(8.0);
% n_rk = 3;
% alfa = [1/2 1/2 1];
% cfl  = 1.95;
% n_rk = 2;
% alfa = [1/2 1];
% cfl  = 1.0;
thresh = -15; % Residual Threshold (log)
%
% Define Nozzle Area
for i = 1:N+2
    A(1,i) = 1.398 + 0.347*tanh(0.8*x_c(1,i)-3.2);
    r(1,i) = sqrt(A(1,i)/pi);
    dAdx(1,i) = 0.8*0.347*sech(0.8*x_c(1,i)-3.2)^2;
end
% figure(1);
% plot(x_c,r,'b');
% hold on;
% plot(x_c,-r,'b');
% axis([0 10 -1.5 1.5]);
% ylabel('r');
% xlabel('x');
% title('Nozzle Profile');
% stop;
%
% File Info
mydir='C:\\Users\\Bucky\\Downloads\\FV2_Results';
subfolder='Ncell'+string(N)+'Jameson';
mkdir(fullfile(mydir,subfolder));
%
% Iterate until convergence
tTot = 0;
niter = 1;
Residual(niter) = 0;
iteration(niter) = 1;
% while Residual(niter) > thresh
for niter = 1:5000
    tStart = tic;
    iteration(niter) = niter;
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
        Qface(1,N+1) = Q(1,N+2);
        Qface(2,N+1) = Q(2,N+2);
        Qface(3,N+1) = Q(3,N+2);
        %
        % Compute Fluxes
        for i = 1:N+1
            F(1,i) = Qface(2,i);
            F(2,i) = (Qface(2,i)^2/Qface(1,i))+((gamma-1)*(Qface(3,i)-0.5*(Qface(2,i)^2/Qface(1,i))));
            F(3,i) = Qface(2,i)*(Qface(3,i)/Qface(1,i) + (gamma-1)*(Qface(3,i)-0.5*Qface(2,i)^2/Qface(1,i))/Qface(1,i));
        end
        %
        % Introduction of Mac-Baldwin Artificial Dissipation
        % for i = 2:N
        %     ubar = 0.5*((abs(u_p1(1,i+1))+c_p1(1,i+1))+(abs(u_p1(1,i))+c_p1(1,i)));
        %     nu = epsilon*ubar*(abs(P_p1(1,i+1)-2*P_p1(1,i)+P_p1(1,i-1))/(P_p1(1,i+1)+2*P_p1(1,i)+P_p1(1,i-1)));
        %     F(:,i) = F(:,i) - nu*(Q(:,i+1)-Q(:,i));
        % end
        %
        % Introduction of Jameson Artificial Dissipation
        for i = 2:N+1
            epsilon2(1,i) = alpha2*(abs(u_p1(1,i))+c_p1(1,i))*(abs(P_p1(1,i+1)-2*P_p1(1,i)+P_p1(1,i-1))/(P_p1(1,i+1)+2*P_p1(1,i)+P_p1(1,i-1)));
        end
        for i = 2:N
            myepsilon2 = (epsilon2(1,i)+epsilon2(1,i+1))/2;
            myepsilon4 = max(0,(alpha4-myepsilon2/(abs(u_p1(1,i))+c_p1(1,i))));
            F(:,i) = F(:,i) - myepsilon2*(Q(:,i+1)-Q(:,i)) + myepsilon4*(Q(:,i+2)-3*Q(:,i+1)+3*Q(:,i)-Q(:,i-1));
        end
        %
        % Compute Source Terms
        S(:,1:N+2) = 0.0;
        for i = 1:N+2 
            S(1,i) = -Q(2,i)*(dAdx(1,i)/A(1,i));
            S(2,i) = -(Q(2,i)^2/Q(1,i))*(dAdx(1,i)/A(1,i));
            S(3,i) = -Q(2,i)*(Q(3,i)/Q(1,i) + ((gamma-1)/Q(1,i))*(Q(3,i)-0.5*Q(2,i)^2/Q(1,i)))*(dAdx(1,i)/A(1,i));
        end
        %
        % Compute Residuals
        for i = 2:N+1
            % R(1,i) = S(1,i) - (F(1,i) - F(1,i-1))/dx;
            % R(2,i) = S(2,i) - (F(2,i) - F(2,i-1))/dx;
            % R(3,i) = S(3,i) - (F(3,i) - F(3,i-1))/dx;
            R(:,i) = S(:,i) - (F(:,i) - F(:,i-1))/dx;
        end
        %
        % Take care of the Residuals at the ends for Q1 and Q_N+2
        F_right(1,1) = Q(2,2);
        F_right(2,1) = (Q(2,2)^2/Q(1,2))*(3-gamma)/2 + (gamma-1)*Q(3,2);
        F_right(3,1) = (gamma*Q(2,2)*Q(3,2)/Q(1,2)) - ((gamma-1)/2)*Q(2,2)^3/Q(1,2)^2;
        R(:,1) = S(:,1) - (F_right(:,1)-F(:,1))/(dx/2);
        %
        F_left(1,1) = Q(2,N+1);
        F_left(2,1) = (Q(2,N+1)^2/Q(1,N+1))*(3-gamma)/2 + (gamma-1)*Q(3,N+1);
        F_left(3,1) = (gamma*Q(2,N+1)*Q(3,N+1)/Q(1,N+1)) - ((gamma-1)/2)*Q(2,N+1)^3/Q(1,N+1)^2;
        R(:,N+2) = S(:,N+2) - (F(:,N+1)-F_left(:,1))/(dx/2);
        %
        % Update Solution
        Q = Q0 + dt*alfa(k)*R;
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
    figure(2);
    plot(x_c,Q(1,:))
    hold on;
    plot(x_c,Q(2,:))
    hold on;
    plot(x_c,Q(3,:))
    legend('Q1', 'Q2', 'Q3')
    axis([0 10 0 6]);
    xlabel('x');
    title('Solved Conservation Variables');
    drawnow;
    hold off;
    % pause(0.01);
    fprintf(1,'Iteration Number = %g Residual = %g\n',niter-1,Residual(niter));
end      

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

figure(1);
plot(x_c,M);
hold on;
plot(x_c,r);
legend('mach number', 'nozzle profile')
axis([0 10 0 2.5]);
ylabel('u');
xlabel('x');
title('Flow Through a Nozzle');

figure(3);
plot(iteration, Residual, 'bo');
grid on;
ylabel('log10(residual)');
xlabel('iteration');
title('Convergence Behavior');

saveas(figure(2),fullfile(mydir,subfolder,subfolder+'_solnCons.jpg'));
saveas(figure(1),fullfile(mydir,subfolder,subfolder+'_soln2.jpg'));
saveas(figure(3),fullfile(mydir,subfolder,subfolder+'_conv.jpg'));

% And Solution Matrices
save(fullfile(mydir,subfolder,subfolder+'_Q.mat'), 'Q')

% And Iteration Info
fid = fopen(fullfile(mydir,subfolder,'iter.txt'),'wt');
fprintf(fid, 'Total Iterations: %s \nTotal CPU-time: %s s', string(niter), string(tTot));
fclose(fid);