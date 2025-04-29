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
% Define variables
clear all; clc;
N = input('Total number of cells = ? ');
x_f = zeros(1,N+1);
x_c = zeros(1,N+1);
Q   = zeros(2,N+2);
F   = zeros(2,N+1);
R   = zeros(2,N+2);
dx = 1.0/N;
x_f(1,1:N+1) = (0:N)/N;
x_c(1,2:N+1) = ((1:N)-0.5)/N;
x_c(1,N+2)   = 1.0;
%n_rk = 4;
%alfa = [1/4 1/3 1/2 1];
%cfl  = sqrt(8.0);
n_rk = 3;
alfa = [1/2 1/2 1];
cfl  = 2.0;
% n_rk = 2;
% alfa = [1/2 1];
% cfl  = 1.0;

dt   = 0.95*dx*cfl;
%
% Iterate until convergence
niter = 1;
Residual(niter) = 0;
while Residual(niter) > -2
    niter = niter + 1;
% for niter = 1:10
    %
    % RK integration
    Q0 = Q;
    for k = 1:n_rk
        %
        % Compute Conservation Variables at Faces
        for i = 2:N
            Qface(1,i) = 0.5*(Q(1,i)+Q(1,i+1));
            Qface(2,i) = 0.5*(Q(2,i)+Q(2,i+1));
        end
        %
        % Compute Fluxes
        for i = 2:N
            F(1,i) = Qface(2,i);
            F(2,i) = Qface(1,i);
        end
        F(1,1) = Q(2,1);
        F(2,1) = Q(1,1);

        F(1,N+1) = Q(2,N+2);
        F(2,N+1) = Q(1,N+2);
        %
        % Compute Source Terms
        S(:,1:N+2) = 0.0;
        for i = 2:N+1
            S(1,i) = 0.0;
            S(2,i) = x_c(1,i);
        end
        % End points
        S(1,1) = 0.0;
        S(2,1) = 1.0;
        %
        % Compute Residuals
        R = 0.0;
        for i = 2:N+1
            R(1,i) = S(1,i) - (F(1,i) - F(1,i-1))/dx;
            R(2,i) = S(2,i) - (F(2,i) - F(2,i-1))/dx;
        end
        %
        % Take care of the Residuals at the ends for Q1 and Q_N+2
        F_Q2_1 = Q(2,2);
        F_Q1_1 = F(1,1);
        F_QN2_1 = F(1,N+1);
        F_QN1_1 = Q(2,N+1);
        F_Q2_2 = Q(1,2);
        F_Q1_2 = F(2,1);
        F_QN2_2 = F(2,N+1);
        F_QN1_2 = Q(1,N+1);
        R(1,1) = S(1,1) - (F_Q2_1-F_Q1_1)/(dx/2);
        R(1,N+2) = S(1,N+2) - (F_QN2_1-F_QN1_1)/(dx/2);
        R(2,1) = S(2,1) - (F_Q2_2-F_Q1_2)/(dx/2);
        R(2,N+2) = S(2,N+2) - (F_QN2_2-F_QN1_2)/(dx/2);
        % R(1,1) = S(1,1) - (F(1,Q2)-F(1,Q1))/(dx/2);
        % R(2,1) = S(2,1) - (F(2,Q2)-F(2,Q1))/(dx/2);
        % R(1,N+2) = S(1,N+2) - (F(1,QN+2)-F(1,QN+1))/(dx/2);
        % R(2,N+2) = S(2,N+2) - (F(2,QN+2)-F(2,QN+1))/(dx/2);
        %
        % Update Solution
        Q = Q0 + dt*alfa(k)*R;
        %
        % Apply BC - Use these for now! I'll go over the derivation
        % 
        % x = 0.0 BCs
        Q(2,1) = Q(2,1) - Q(1,1);
        Q(1,1) = 0.0;
        %
        % x = 1.0 BCs
        Q(1,N+2) = 0.5*(Q(1,N+2) + Q(2,N+2));
        Q(2,N+2) = Q(1,N+2);
    end
%
    Residual(niter) = log10(norm(R)/(2*(N+2)));           
    figure(1);
    plot(x_c,Q);
    axis([0 1 0 0.6]);
    fprintf(1,'Iteration Number = %g Residual = %g\n',niter,Residual(niter));
end