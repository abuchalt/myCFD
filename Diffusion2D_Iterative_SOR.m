%% Diffusion2D_Iterative
% ------------------------------------------------------------------------------
% This is an iterative finite-difference solver for 2-dimensional, Laplacian 
% heat diffusion in cartesian coordinates using Dirichlet boundary conditions
% ∇⋅∇ T = 0
% This solver implements a successive over-relaxation algorithm by introducing a
% pseudo-timestepping term in order to converge more quickly upon the solution
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

% Define temperature at boundry
T_boundry = 100.0;

% Define iteration parameters
residual = 1.0E5;
epsilon = 1.0E-9;

% Define relaxation factor (pseudo-timestepping)
Deltatau = 1.8;

% Calculate step sizes
Deltax = w/(i_max-1);
Deltay = h/(j_max-1);

% Define x and y values
for i = 1:i_max
    for j = 1:j_max
        x(i,j) = Deltax*(i-1);
        y(i,j) = Deltay*(j-1);
        % k = pmap(i,j,i_max); % Alternative to reshape (below)
        % T_plot(i,j) = T(k,1);
    end
end

%% Script
% ------------------------------------------------------------------------------

% Discretize and convert to a linear system A*T = b

% Init
T = zeros(i_max,j_max);
% Apply BCs
T(1,:) = T_boundry;

% Define Coefficient Matrix
iter = 1;
T_old = T;
tstart = tic; % Measure computation time
while (residual > epsilon)  
    for i = 2:i_max-1
        for j = 2:j_max-1
            % T(i,j) = ((T(i+1,j)+T(i-1,j))/Deltax^2 + (T(i,j+1)+T(i,j-1))/Deltay^2)/(2.0/Deltax^2 + 2.0/Deltay^2); % Unrelaxed Gauss-Seidal
            % T(i,j) = T(i,j) + ( (T(i+1,j)-2*T(i,j)+T(i-1,j))/Deltax^2 + (T(i,j+1)-2*T(i,j)+T(i,j-1))/Deltay^2 ) / (2.0/Deltax^2 + 2.0/Deltay^2); % Relaxation form w/ same performance (no OVER-relaxation)
            % overrelaxation: Δτ > 1.0   --Increases rate of convergence, too large will diverge
            T(i,j) = T(i,j) + Deltatau*( (T(i+1,j)-2*T(i,j)+T(i-1,j))/Deltax^2 + (T(i,j+1)-2*T(i,j)+T(i,j-1))/Deltay^2 ) / (2.0/Deltax^2 + 2.0/Deltay^2);
        end
    end
    iter = iter + 1;
    residual = norm(T-T_old);
    fprintf(1,'iter = %i, residual = %g\n',iter,log10(residual))
    T_old = T;
    % draw pseudo-time evolution live
    % contour(x,y,T,'LineWidth',2);
    % drawnow;
end
toc(tstart);

%% Output Results
% ------------------------------------------------------------------------------

% Produce isotherm level contour plot of data
contour(x,y,T,'LineWidth',2)