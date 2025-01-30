%% Diffusion2D
% ------------------------------------------------------------------------------
% This is a finite-difference solver for 2-dimensional, Laplacian heat diffusion
% in cartesian coordinates using Dirichlet boundary conditions
% ∇⋅∇ T = 0
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

% Calculate step sizes
Deltax = w/(i_max-1);
Deltay = h/(j_max-1);

%% Script
% ------------------------------------------------------------------------------

% Discretize and convert to a linear system A*T = b

% Init
% A = eye(i_max*j_max); % Coefficient Matrix
A = spalloc(i_max*j_max, 5*i_max, j_max); % Allocate Coeff Matrix Sparsely
for i = 1:i_max*j_max
    A(i,i) = 1.0;
end
b = zeros(i_max*j_max,1);

% Define Coefficient Matrix
for i = 2:i_max-1
    for j = 2:j_max-1
        k = pmap(i, j, i_max);
        k_e = k + 1;
        k_w = k - 1;
        k_n = k + i_max;
        k_s = k - i_max;

        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        A(k,k) = -2.0/Deltax^2 - 2.0/Deltay^2;
        A(k,k_e) = 1.0/Deltax^2;
        A(k,k_w) = 1.0/Deltax^2;
        A(k,k_n) = 1.0/Deltay^2;
        A(k,k_s) = 1.0/Deltay^2;
    end
end
A = sparse(A); % Define A sparse

% Apply BCs
for i = 1:1
    for j = 1:j_max
        k = pmap(i,j,i_max);
        b(k,1) = T_boundry;
    end
end

% Solve the linear system via mldivide
T = A\b;

%% Output Results
% ------------------------------------------------------------------------------

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
T_plot = reshape(T,i_max,j_max);

% Produce isotherm level contour plot of data
contour(x,y,T_plot)

%% Functions
% ------------------------------------------------------------------------------

function k = pmap(i, j, i_max)
    k = i + (j-1)*i_max;
end

function [i, j] = revPmap(k, i_max)
    i = 1 + mod((k-1), i_max);
    j = 1 + ((k-i)/i_max);
end