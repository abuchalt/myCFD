%% newtonsMethod
% ------------------------------------------------------------------------------
% This is a newton's method solver for pairs of nonlinear functions of
% 2 independent variables
% ------------------------------------------------------------------------------
clear all; close all; clc;

% x^2 + y^2 = 25
% 3xy + 2y^2 = 68

% Initial Guess
x_0 = 1;
y_0 = 2;

% Define Residual + Convergence Criterion
residual = 1.0E5; % init residual
epsilon = 1.0E-12; % drive residual down to this value before terminating

% Initialize iteration counter
iter = 0;

% Use a while loop
while(residual>epsilon)
    % Build the Jacobian
    J = [2*x_0, 2*y_0;
         3*y_0, 3*x_0+4*y_0];
    % And the RHS
    F = [25-x_0^2-y_0^2;
         68-3*x_0*y_0-2*y_0^2];
    % Solve for deltas
    DeltaX = J\F;
    % Determine new `guess'
    x_0 = x_0 + DeltaX(1);
    y_0 = y_0 + DeltaX(2);
    % Calculate residual from deltas
    residual = norm(DeltaX)/2;
    iter = iter + 1;
    fprintf(1,'iter = %i, log10(residual) = %g\n',iter,log10(residual));
end

fprintf(1,'x = %g\n', x_0);
fprintf(1,'y = %g\n', y_0);