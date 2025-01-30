% Computational Parameters
% ------------------------------------------------------------------------------

i_max = 5;
j_max = 5;

% Script
% ------------------------------------------------------------------------------

% Solving the Laplacian Heat Diffusion Equation: ∇⋅∇ T = 0

% Discretize and convert to a linear system A*T = B

% Solve the system via mldivide
% T = A\B

% Functions
% ------------------------------------------------------------------------------

function k = seq(i, j, i_max)
    k = i + (j-1)*i_max;
end

function [i, j] = revSeq(k, i_max)
    i = 1 + mod((k-1), i_max);
    j = 1 + ((k-i)/i_max);
end