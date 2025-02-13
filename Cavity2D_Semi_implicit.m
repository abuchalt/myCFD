%% Cavity2D_Semi_implicit
% ------------------------------------------------------------------------------
% This is a semi-implicit solver for 2-dimensional Navier-Stokes equations, by
% lagging the velocity terms in the vorticity transport equation, for a cavity
% with tangential flow
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% Computational Parameters
% ------------------------------------------------------------------------------

% Define physical domain
h = 1.0;
w = 1.0;

% Define mesh size
fprintf('Maximum number of points in x-direction') % separate print and input
                                                   % b/c vscode extension
i_max = input('');
fprintf('Maximum number of points in y-direction')
j_max = input('');

% Input parameters 
Re = 100.0; % Reynold's number (kinematic viscosity)
u_lid = 1.0; % velocity at top boundry

% Calculate step sizes
Deltax = w/(i_max-1);
Deltay = h/(j_max-1);

% Define pseudo-timestep
% Equivalent to parabolic equation in 1D - 
% the diffusion number determined by Reynold's number for this scheme should be s.t. d <= 0.5 to ensure stability 
% d = 1/Re * (dtau/dx)^2
Deltatau = Re*Deltax^2/2.0;
% And define variables for pseudo-timestepping
residual = 1.0E5; % init residual
epsilon = 1.0E-12; % drive residual down to this value before terminating

%% Script
% ------------------------------------------------------------------------------

% Discretize and convert to a linear system A*T = b

% Init
A_Psi = spalloc(i_max*j_max, 5*i_max, j_max); % Allocate Streamfxn Coeff Matrix Sparsely
A_Omega = spalloc(i_max*j_max, 5*i_max, j_max); % Allocate Vorticity Coeff Matrix Sparsely
for i = 1:i_max*j_max
    A_Psi(i,i) = 1.0;
    A_Omega(i,i) = 1.0;
end
b_Psi = zeros(i_max*j_max,1); % init RHS
b_Omega = zeros(i_max*j_max,1);

% Define Solution Variables (1D because we use pointer mapping)
Psi = zeros(i_max*j_max,1);
Omega = zeros(i_max*j_max,1);
u = zeros(i_max*j_max,1);
v = zeros(i_max*j_max,1);
% "Old" Solution for finding residual
Psi_old = zeros(i_max*j_max,1);
Omega_old = zeros(i_max*j_max,1);

% Define Coefficient Matrix
for i = 2:i_max-1
    for j = 2:j_max-1
        k = pmap(i, j, i_max);
        k_e = k + 1;
        k_w = k - 1;
        k_n = k + i_max;
        k_s = k - i_max;

        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        A_Psi(k,k) = 1.0/Deltatau + (2.0/Deltax^2 + 2.0/Deltay^2)/Re;
        A_Psi(k,k_e) = (-1.0/Deltax^2)/Re;
        A_Psi(k,k_w) = (-1.0/Deltax^2)/Re;
        A_Psi(k,k_n) = (-1.0/Deltay^2)/Re;
        A_Psi(k,k_s) = (-1.0/Deltay^2)/Re;

        % assemble RHS
        b_Psi(k,1) = Psi(k,1)/Deltatau + Omega(k,1)/Re;
    end
end
A_Psi = sparse(A_Psi); % Enforce A_Psi sparse

% Apply 

% Begin iteration
iter = 0;
while (residual > epsilon)

    % Enforce velocity boundry condition
    j = j_max
    for i = 1:i_max
        k = pmap(i, j, i_max);
        u(k,1) = u_lid;
    end

    for i = 2:i_max-1
        for j = 2:j_max-1
            k = pmap(i, j, i_max);
            k_e = k + 1;
            k_w = k - 1;
            k_n = k + i_max;
            k_s = k - i_max;

            % Define velocity (u,v) based on initial Psi values
            u(k,1) = (Psi(k_n,1)-Psi(k_s,1))/(2.0*Deltay);
            v(k,1) = -(Psi(k_w,1)-Psi(k_e,1))/(2.0*Deltax);

            % Update vorticity coefficient matrix
            A_Omega(k,k) = 1.0/Deltatau + (2.0/Deltax^2 + 2.0/Deltay^2)/Re;
            A_Omega(k,k_e) = (-1.0/Deltax^2)/Re + u(k,1)/(2.0*Deltax);
            A_Omega(k,k_w) = (-1.0/Deltax^2)/Re - u(k,1)/(2.0*Deltax);
            A_Omega(k,k_n) = (-1.0/Deltay^2)/Re + v(k,1)/(2.0*Deltay);
            A_Omega(k,k_s) = (-1.0/Deltay^2)/Re - v(k,1)/(2.0*Deltay);

            b_Psi(k,1) = Psi(k,1)/Deltatau + Omega(k,1)/Re;
            b_Omega(k,1) = Omega(k,1)/Deltatau;
        end
    end

    % Apply BCs
    % [INSERT CODE HERE]
    % for i = 1:1
    %     for j = 1:j_max
    %         k = pmap(i,j,i_max);
    %         b(k,1) = T_boundry;
    %     end
    % end
    
    % Solve the linear systems for Psi and Omega
    Psi = A_Psi\b_Psi;
    Omega = A_Omega\b_Omega;

    % Compute the new residual
    residual = 0.0;
    for i = 1:i_max
        for j = 1:j_max
            k = pmap(i, j, i_max);
            residual = residual + (Psi(k,1) - Psi_old(k,1))^2 + (Omega(k,1) - Omega_old(k,1))^2;
        end
    end
    residual = sqrt(residual)/2.0; % Normalize for # of equations being solved
    residual = residual/(i_max*j_max); % Normalize for DOF
    iter = iter+1;
    Psi_old = Psi; % Update "old" values
    Omega_old = Omega;
end

%% Output Results
% ------------------------------------------------------------------------------

% [FIX THIS TOO]

% Init

% Define x and y values
for i = 1:i_max
    for j = 1:j_max
        x(i,j) = Deltax*(i-1);
        y(i,j) = Deltay*(j-1);
        % k = pmap(i,j,i_max); % Alternative to reshape (below)
        % T_plot(i,j) = T(k,1);
    end
end

% Reshape Solution T to imax/jmax mesh
% T_plot = reshape(T,i_max,j_max);

% Produce isotherm level contour plot of data
% contour(x,y,T_plot)

%% Functions
% ------------------------------------------------------------------------------

function k = pmap(i, j, i_max)
    k = i + (j-1)*i_max;
end

% function [i, j] = revPmap(k, i_max)
%     i = 1 + mod((k-1), i_max);
%     j = 1 + ((k-i)/i_max);
% end