%% Cavity2D_Semi_implicit_Fully_coupled
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
% the diffusion number determined by Reynold's number for an UNCOUPLED scheme should be s.t. d <= 0.5 to ensure stability 
% d = 1/Re * (dtau/dx)^2
scale = 0.05;
Deltatau = Re*Deltax^2/scale;
% And define variables for pseudo-timestepping
residual = 1.0E5; % init residual
epsilon = 1.0E-12; % drive residual down to this value before terminating

% Define x and y values in spatial domain
for i = 1:i_max
    for j = 1:j_max
        x(i,j) = Deltax*(i-1);
        y(i,j) = Deltay*(j-1);
    end
end

% File Info
mydir='C:\\Users\\Bucky\\Downloads\\2DCavity_Results';
subfolder='Re'+string(Re)+'_'+string(i_max)+'x'+string(j_max);
mkdir(fullfile(mydir,subfolder));

%% Script
% ------------------------------------------------------------------------------

% Discretize and convert to a linear system A*T = b

% Init
A_PsiPsi = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Allocate Streamfxn Coeff Matrix Sparsely
A_PsiOmega = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Allocate Coeff Matrix for Streamfxn Dependence on Vorticity
A_OmegaOmega = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Allocate Vorticity Coeff Matrix Sparsely
A_OmegaPsi = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Allocate Coeff Matrix for Vorticity Dependence on Streamfxn
for i = 1:i_max*j_max
    A_PsiPsi(i,i) = 1.0;
    A_OmegaOmega(i,i) = 1.0;
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
        A_PsiPsi(k,k) = (1.0/Deltatau) + (2.0/Deltax^2)/Re + (2.0/Deltay^2)/Re;
        A_PsiPsi(k,k_e) = (-1.0/Deltax^2)/Re;
        A_PsiPsi(k,k_w) = (-1.0/Deltax^2)/Re;
        A_PsiPsi(k,k_n) = (-1.0/Deltay^2)/Re;
        A_PsiPsi(k,k_s) = (-1.0/Deltay^2)/Re;

        A_PsiOmega(k,k) = -1.0/Re;
    end
end
A_PsiPsi = sparse(A_PsiPsi); % Enforce Coeffs sparse
A_PsiOmega = sparse(A_PsiOmega); 
A_OmegaOmega = sparse(A_OmegaOmega); 
A_OmegaPsi = sparse(A_OmegaPsi);

% Begin iteration
tTot = 0;
iter = 0;
while (residual > epsilon)

    tStart = tic;

    % Enforce velocity boundry condition
    j = j_max;
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
            v(k,1) = -(Psi(k_e,1)-Psi(k_w,1))/(2.0*Deltax);

            % Update vorticity coefficient matrix
            A_OmegaOmega(k,k) = 1.0/Deltatau + (2.0/Deltax^2 + 2.0/Deltay^2)/Re;
            A_OmegaOmega(k,k_e) = (-1.0/Deltax^2)/Re + u(k,1)/(2.0*Deltax);
            A_OmegaOmega(k,k_w) = (-1.0/Deltax^2)/Re - u(k,1)/(2.0*Deltax);
            A_OmegaOmega(k,k_n) = (-1.0/Deltay^2)/Re + v(k,1)/(2.0*Deltay);
            A_OmegaOmega(k,k_s) = (-1.0/Deltay^2)/Re - v(k,1)/(2.0*Deltay);

            b_Psi(k,1) = Psi(k,1)/Deltatau;
            b_Omega(k,1) = Omega(k,1)/Deltatau;
        end
    end

    % Apply BCs
    % Left BC
    for i = 1:1
        for j = 1:j_max
            k = pmap(i,j,i_max);
            k_e = k + 1;
            k_ee = k + 2;
            % Vorticity via second-order forward difference
            A_OmegaOmega(k,k) = -1.0/Re;
            A_OmegaPsi(k,k_e) = (-8.0/(2.0*(Deltay^2)))/Re;
            A_OmegaPsi(k,k_ee) = (1.0/(2.0*(Deltay^2)))/Re;
            b_Omega(k,1) = 0.0;
            % No-slip (streamfxn = 0)
            A_PsiPsi(k,k) = 1.0;
            b_Psi(k,1) = 0.0;
        end
    end
    % Right BC
    for i = i_max:i_max
        for j = 1:j_max
            k = pmap(i,j,i_max);
            k_w = k - 1;
            k_ww = k - 2;
            % Vorticity via second-order backward difference
            A_OmegaOmega(k,k) = -1.0/Re;
            A_OmegaPsi(k,k_w) = (-8.0/(2.0*(Deltay^2)))/Re;
            A_OmegaPsi(k,k_ww) = (1.0/(2.0*(Deltay^2)))/Re;
            b_Omega(k,1) = 0.0;
            % No-slip (streamfxn = 0)
            A_PsiPsi(k,k) = 1.0;
            b_Psi(k,1) = 0.0;
        end
    end
    % Bottom BC
    for i = 1:i_max
        for j = 1:1
            k = pmap(i,j,i_max);
            k_n = k + i_max;
            k_nn = k + 2*i_max;
            % Vorticity via second-order forward difference
            A_OmegaOmega(k,k) = -1.0/Re;
            A_OmegaPsi(k,k_n) = (-8.0/(2.0*(Deltay^2)))/Re;
            A_OmegaPsi(k,k_nn) = (1.0/(2.0*(Deltay^2)))/Re;
            b_Omega(k,1) = 0.0;
            % No-slip (streamfxn = 0)
            A_PsiPsi(k,k) = 1.0;
            b_Psi(k,1) = 0.0;
        end
    end
    % Top BC
    for i = 1:i_max
        for j = j_max:j_max
            k = pmap(i,j,i_max);
            k_s = k - i_max;
            k_ss = k - 2*i_max;
            % Vorticity via second-order backward difference with lid-driven BC
            A_OmegaOmega(k,k) = -1.0/Re;
            A_OmegaPsi(k,k_s) = (-8.0/(2.0*(Deltay^2)))/Re;
            A_OmegaPsi(k,k_ss) = (1.0/(2.0*(Deltay^2)))/Re;
            b_Omega(k,1) = (3.0*u_lid/Deltay)/Re;
            % Lid Velocity enforced (streamfxn = constant)
            A_PsiPsi(k,k) = 1.0;
            b_Psi(k,1) = 0.0;
        end
    end
    
    % Solve the linear systems for Psi and Omega
    M = [A_PsiPsi, A_PsiOmega;
         A_OmegaPsi, A_OmegaOmega];
    b = [b_Psi;
         b_Omega];
    M = sparse(M); % Enforce Sparse
    sol = M\b;
    Psi = sol(1:i_max*j_max);
    Omega = sol((i_max*j_max)+1:2*i_max*j_max);

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
    Psi_old = Psi; % Update "old" values
    Omega_old = Omega;

    tTot = tTot + toc(tStart);

    % Plot solution
    if mod(iter,100) == 0
        uplot = reshape(u, i_max, j_max);
        vplot = reshape(v, i_max, j_max);
        figure(1);
        % Plot level curves for vorticity - set levels according to Ghia et al.
        contour(x,y,reshape(Omega, i_max, j_max),[-5.0 -4.0 -3.0 -2.0 -1.0 -0.5 0.0 0.5 1.0 2.0 3.0],'LineWidth',2.0);
        ylabel('y');
        xlabel('x');
        title('Vorticity Contour');
        pbaspect([w h 1]);
        figure(2);
        % Plot level curve for streamfxn - set levels according to Ghia et al.
        contour(x,y,reshape(Psi, i_max, j_max),[-0.1175 -0.1150 -0.11 -0.1 -0.09 -0.07 -0.05 -0.03 -0.01 -1E-4 -1E-5 -1E-7 -1E-10 1E-8 1E-7 1E-6 1E-5 5E-5 1E-4 2.5E-4 5E-4 1E-3 1.5E-3 3E-3], 'LineWidth',2.0)
        ylabel('y');
        xlabel('x');
        title('Streamline Pattern');
        pbaspect([w h 1]);
        figure(3);
        % Plot u-v vector field (velocity)
        quiver(x,y,uplot,vplot,20);
        ylabel('y');
        xlabel('x');
        title('Velocity Vector Field');
        axis([0 1 0 1]);
        pbaspect([w h 1]);
        figure(4);
        hold on;
        % Plot log10(residual)
        plot(iter,log10(residual),'bo');
        grid on;
        ylabel('log10(residual)');
        xlabel('iteration');
        title('Convergence Behavior');
        hold off;
        drawnow;
    end
    fprintf(1,'iter = %i, residual = %g\n',iter,log10(residual));
    iter = iter+1;
end

%% Output Results
% ------------------------------------------------------------------------------

figure(5);
plot(x(:,round((j_max+1)/2)), vplot(:,round((j_max+1)/2)));
grid on;
ylabel('v');
xlabel('x');
title('Vertical Component of Velocity Through Geometric Center');

figure(6);
plot(uplot(round((i_max+1)/2),:), y(round((i_max+1)/2),:));
grid on;
ylabel('y');
xlabel('u');
title('Horizontal Component of Velocity Through Geometric Center');

% Save Final Plots
saveas(figure(1),fullfile(mydir,subfolder,subfolder+'_vorticity.jpg'));
saveas(figure(2),fullfile(mydir,subfolder,subfolder+'_streamfxn.jpg'));
saveas(figure(3),fullfile(mydir,subfolder,subfolder+'_velocity.jpg'));
saveas(figure(4),fullfile(mydir,subfolder,subfolder+'_convergence.jpg'));
saveas(figure(5),fullfile(mydir,subfolder,subfolder+'_verticalVel.jpg'));
saveas(figure(6),fullfile(mydir,subfolder,subfolder+'_horizontalVel.jpg'));

% And Solution Matrices
save(fullfile(mydir,subfolder,subfolder+'_Psi.mat'), 'Psi')
save(fullfile(mydir,subfolder,subfolder+'_Omega.mat'), 'Omega')

% And Timing Info
tAvg = tTot/iter;
fid = fopen(fullfile(mydir,subfolder,'time.txt'),'wt');
fprintf(fid, 'Total CPU-time: %s s\nAverage Time per Iteration: %s s\nPseudotime Step: Re*Deltax^2/%s', string(tTot), string(tAvg), string(CFL));
fclose(fid);

%% Functions
% ------------------------------------------------------------------------------

function k = pmap(i, j, i_max)
    k = i + (j-1)*i_max;
end