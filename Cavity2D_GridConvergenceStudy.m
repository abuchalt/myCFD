%% Cavity2D_GridConvergenceStudy
% ------------------------------------------------------------------------------
% This is a script for estimating and visualizing convergence of truncation
% error for the lid-driven cavity problem with square, uniform meshes
% ------------------------------------------------------------------------------
clear all; close all; clc;

%% User Input
% ------------------------------------------------------------------------------

% Params
Re = 100;
meshSizes = [257, 129, 65];

% File Info
mydir='C:\\Users\\Bucky\\Downloads\\2DCavity_Results';

%% Import Data
% ------------------------------------------------------------------------------

% Quick Maths
meshSizes = sort(meshSizes);
minN = meshSizes(1);
meshes = size(meshSizes, 2);
assert(meshes>=3, 'Need at Least 3 Grids')
midline = (minN+1)/2;

% Init
OmegaPages = zeros(minN,minN,meshes);
PsiPages = zeros(minN,minN,meshes);

% Evaluate Error only at coincident points of uniform meshes
% function E = getError(f_fine, f_coarse)
%     Nfine = (size(f_fine,1)-1)
%     Ncoarse = (size(f_coarse,1)-1)
%     r = round(Nfine/Ncoarse);
%     f_fineresize = f_fine(1:r:Nfine);
%     E = size(f_fineresize,1);
% end

% Store Resized Solution Matrices
for i = 1:meshes
    N = meshSizes(i);
    r = (N-1)/(minN-1);
    assert(floor(r)==r, 'Grid Refinement Ratio is Non-Integer');
    subfolder='Re'+string(Re)+'_'+string(N)+'x'+string(N);
    myOmega = matfile(fullfile(mydir,subfolder,subfolder+'_Omega.mat'));
    myPsi = matfile(fullfile(mydir,subfolder,subfolder+'_Psi.mat'));
    Omega = reshape(myOmega.Omega, N, N);
    OmegaPages(1:minN,1:minN,i) = Omega(1:r:N,1:r:N);
    Psi = reshape(myPsi.Psi, N, N);
    PsiPages(1:minN,1:minN,i) = Psi(1:r:N,1:r:N);
end

%% Estimate Order of Convergence by Vorticity at Center
% ------------------------------------------------------------------------------

% Init
f = zeros(meshes, 1);

% Evaluate Order of Convergence Directly
for i = 1:meshes
    f(i, 1) = OmegaPages(midline,midline,i);
end
r = 2;
p1 = log(norm(f(1,1)-f(2,1))/norm(f(2,1)-f(3,1)))/log(r);
fprintf('Order of Grid Convergence by Vorticity at Center: %g\n', p1);

% Demonstrate Order of Convergence Graphically
h = 1./meshSizes;
% Transform Into Expected Linear Behavior Space
h2 = h.^2;
magOmega = abs(f(:,1));
% Line of Fit
coefficients = polyfit(h2, magOmega, 1);
xFit = linspace(min(h2), max(h2), 3);
yFit = polyval(coefficients , xFit);
% Inspect Closeness
figure(1);
plot(h2, magOmega, 'o'); % Plot Real Data
hold on;
plot(xFit, yFit, '--'); % Plot Fit
hold off;
ylabel('Magnitude of Vorticity at Mesh Center','interpreter','latex');
xlabel('Square of Relative Grid Spacing $h^2$','interpreter','latex');
title('Estimation of Order of Convergence by Vorticity for Re='+string(Re),'interpreter','latex');
filename = 'GridConvergence_Vorticity_Re'+string(Re)+'.jpg';
saveas(figure(1),fullfile(mydir,filename));

%% Estimate Order of Convergence by Streamfunction at Center
% ------------------------------------------------------------------------------

% Init
f = zeros(meshes, 1);

% Evaluate Order of Convergence Directly
for i = 1:meshes
    f(i, 1) = PsiPages(midline,midline,i);
end
r = 2;
p2 = log(norm(f(1,1)-f(2,1))/norm(f(2,1)-f(3,1)))/log(r);
fprintf('Order of Grid Convergence by Streamfxn at Center: %g\n', p2);

% Demonstrate Order of Convergence Graphically
h = 1./meshSizes;
% Transform Into Expected Linear Behavior Space
h2 = h.^2;
magOmega = abs(f(:,1));
% Line of Fit
coefficients = polyfit(h2, magOmega, 1);
xFit = linspace(min(h2), max(h2), 3);
yFit = polyval(coefficients , xFit);
% Inspect Closeness
figure(2);
plot(h2, magOmega, 'o'); % Plot Real Data
hold on;
plot(xFit, yFit, '--'); % Plot Fit
hold off;
ylabel('Magnitude of Streamfunction at Mesh Center','interpreter','latex');
xlabel('Square of Relative Grid Spacing $h^2$','interpreter','latex');
title('Estimation of Order of Convergence by Streamfunction for Re='+string(Re),'interpreter','latex');
filename = 'GridConvergence_Streamfxn_Re'+string(Re)+'.jpg';
saveas(figure(2),fullfile(mydir,filename));

%% Estimate Order of Convergence by Global Relative Error
% ------------------------------------------------------------------------------

E = zeros(meshes, 1);
for i = 1:meshes-1
    E(i, 1) = (norm(PsiPages(:,:,i) - PsiPages(:,:,meshes))+norm(OmegaPages(:,:,i) - OmegaPages(:,:,meshes)))/(minN^2*2);
end

% Estimate Order of Convergence Graphically
h = 1./meshSizes;
% Line of Fit
fit = fit(reshape(h, meshes, 1),E,'C*x^p');
fprintf('Order of Grid Convergence by Global Relative Error: %g\n', fit.p);
% Inspect Closeness
figure(3);
plot(h, E, 'o'); % Plot Real Data
hold on;
plot(fit, '--'); % Plot Fit
hold off;
ylabel('Magnitude of Relative Global Error','interpreter','latex');
xlabel('Relative Grid Spacing $h$','interpreter','latex');
title('Estimation of Order of Convergence by Global Relative Error for Re='+string(Re),'interpreter','latex');
filename = 'GridConvergence_GlobalError_Re'+string(Re)+'.jpg';
saveas(figure(3),fullfile(mydir,filename));

%% Save Values
% ------------------------------------------------------------------------------
filename = 'GridConvergence_Re'+string(Re)+'.txt';
fid = fopen(fullfile(mydir,filename),'wt');
fprintf(fid, 'Order of Grid Convergence by Vorticity at Center: %g\nOrder of Grid Convergence by Streamfxn at Center: %g\nOrder of Grid Convergence by Global Relative Error: %g\n', p1, p2, fit.p);
fclose(fid);