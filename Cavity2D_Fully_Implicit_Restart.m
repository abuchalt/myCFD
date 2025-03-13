%% Cavity2D_Fully_Implicit_Restart
% ------------------------------------------------------------------------------
% This is a fully implicit solver for 2-dimensional Navier-Stokes equations, via
% perturbation method applied to the vorticity transport equation, for a cavity
% with tangential flow
%
% This script allows the iteration to be paused and continued from a new 
% pseudotimestep in order to aid the approach to high Reynolds Number solutions
% ------------------------------------------------------------------------------

% Re-init residual to restart iteration
residual = 1.0E5;

% Define Coefficient Matrix
for i = 2:i_max-1
    for j = 2:j_max-1
        k = pmap(i, j, i_max);
        k_e = k + 1;
        k_w = k - 1;
        k_n = k + i_max;
        k_s = k - i_max;

        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        A_PsiPsi(k,k) = (1.0/Deltatau) + (2.0/Deltax^2 + 2.0/Deltay^2)/Re;
        % A_PsiPsi(k,k) = (2.0/Deltax^2 + 2.0/Deltay^2)/Re;
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

            % Define velocity (u,v) based on current Psi values
            u(k,1) = (Psi(k_n,1)-Psi(k_s,1))/(2.0*Deltay);
            v(k,1) = -(Psi(k_e,1)-Psi(k_w,1))/(2.0*Deltax);

            % Update vorticity coefficient matrix
            A_OmegaOmega(k,k) = (1.0/Deltatau) + (2.0/Deltax^2 + 2.0/Deltay^2)/Re;
            % A_OmegaOmega(k,k) = (2.0/Deltax^2 + 2.0/Deltay^2)/Re;
            A_OmegaOmega(k,k_e) = (-1.0/Deltax^2)/Re + u(k,1)/(2.0*Deltax);
            A_OmegaOmega(k,k_w) = (-1.0/Deltax^2)/Re - u(k,1)/(2.0*Deltax);
            A_OmegaOmega(k,k_n) = (-1.0/Deltay^2)/Re + v(k,1)/(2.0*Deltay);
            A_OmegaOmega(k,k_s) = (-1.0/Deltay^2)/Re - v(k,1)/(2.0*Deltay);

            A_OmegaPsi(k,k_e) = -(Omega(k_n,1)-Omega(k_s,1))/(4.0*Deltay*Deltax);
            A_OmegaPsi(k,k_w) = (Omega(k_n,1)-Omega(k_s,1))/(4.0*Deltay*Deltax);
            A_OmegaPsi(k,k_n) = (Omega(k_e,1)-Omega(k_w,1))/(4.0*Deltay*Deltax);
            A_OmegaPsi(k,k_s) = -(Omega(k_e,1)-Omega(k_w,1))/(4.0*Deltay*Deltax);

            b_Omega(k,1) = -((Psi(k_n,1)-Psi(k_s,1))/(2.0*Deltay))*((Omega(k_e,1)-Omega(k_w,1))/(2.0*Deltax)) + ((Psi(k_e,1)-Psi(k_w,1))/(2.0*Deltax))*((Omega(k_n,1)-Omega(k_s,1))/(2.0*Deltay)) + (((Omega(k_e,1)-2*Omega(k,1)+Omega(k_w,1))/(Deltax^2)) + ((Omega(k_n,1)-2*Omega(k,1)+Omega(k_s,1))/(Deltay^2)))/Re;
            % b_Omega(k,1) = (-(Psi(k_n,1)*Omega(k_e,1)) + (Psi(k_n,1)*Omega(k_w,1)) + (Psi(k_s,1)*Omega(k_e,1)) - (Psi(k_s,1)*Omega(k_w,1)) + (Omega(k_n,1)*Psi(k_e,1)) - (Omega(k_n,1)*Psi(k_w,1)) - (Omega(k_s,1)*Psi(k_e,1)) + (Omega(k_s,1)*Psi(k_w,1)))/(4.0*Deltax*Deltay) + (((Omega(k_e,1)-2*Omega(k,1)+Omega(k_w,1))/(Deltax^2)) + ((Omega(k_n,1)-2*Omega(k,1)+Omega(k_s,1))/(Deltay^2)))/Re;
            b_Psi(k,1) = (Omega(k,1) + (Psi(k_e,1)-2*Psi(k,1)+Psi(k_w,1))/(Deltax^2) + (Psi(k_n,1)-2*Psi(k,1)+Psi(k_s,1))/(Deltay^2))/Re;
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
            A_OmegaOmega(k,k) = -1.0;
            A_OmegaPsi(k,k_e) = (-8.0/(2.0*(Deltax^2)));
            A_OmegaPsi(k,k_ee) = (1.0/(2.0*(Deltax^2)));
            b_Omega(k,1) = (Omega(k,1) + (8*Psi(k_e)-Psi(k_ee))/(2.0*(Deltax^2)));
            % No-slip (streamfxn = constant 0)
            A_PsiPsi(k,k) = 1.0;
            b_Psi(k,1) = -Psi(k,1);
        end
    end
    % Right BC
    for i = i_max:i_max
        for j = 1:j_max
            k = pmap(i,j,i_max);
            k_w = k - 1;
            k_ww = k - 2;
            % Vorticity via second-order backward difference
            A_OmegaOmega(k,k) = -1.0;
            A_OmegaPsi(k,k_w) = (-8.0/(2.0*(Deltax^2)));
            A_OmegaPsi(k,k_ww) = (1.0/(2.0*(Deltax^2)));
            b_Omega(k,1) = (Omega(k,1) + (8*Psi(k_w)-Psi(k_ww))/(2.0*(Deltax^2)));
            % No-slip (streamfxn = constant 0)
            A_PsiPsi(k,k) = 1.0;
            b_Psi(k,1) = -Psi(k,1);
        end
    end
    % Bottom BC
    for i = 1:i_max
        for j = 1:1
            k = pmap(i,j,i_max);
            k_n = k + i_max;
            k_nn = k + 2*i_max;
            % Vorticity via second-order forward difference
            A_OmegaOmega(k,k) = -1.0;
            A_OmegaPsi(k,k_n) = (-8.0/(2.0*(Deltay^2)));
            A_OmegaPsi(k,k_nn) = (1.0/(2.0*(Deltay^2)));
            b_Omega(k,1) = (Omega(k,1) + (8*Psi(k_n)-Psi(k_nn))/(2.0*(Deltay^2)));
            % No-slip (streamfxn = constant 0)
            A_PsiPsi(k,k) = 1.0;
            b_Psi(k,1) = -Psi(k,1);
        end
    end
    % Top BC
    for i = 1:i_max
        for j = j_max:j_max
            k = pmap(i,j,i_max);
            k_s = k - i_max;
            k_ss = k - 2*i_max;
            % Vorticity via second-order backward difference with lid-driven BC
            A_OmegaOmega(k,k) = -1.0;
            A_OmegaPsi(k,k_s) = (-8.0/(2.0*(Deltay^2)));
            A_OmegaPsi(k,k_ss) = (1.0/(2.0*(Deltay^2)));
            b_Omega(k,1) = (Omega(k,1) + (3.0*u_lid/Deltay) + (8*Psi(k_s)-Psi(k_ss))/(2.0*(Deltay^2)));
            % Lid Velocity enforced (streamfxn = constant)
            A_PsiPsi(k,k) = 1.0;
            b_Psi(k,1) = -Psi(k,1);
        end
    end
    
    % Solve the linear systems for Psi and Omega
    M = [A_PsiPsi, A_PsiOmega;
         A_OmegaPsi, A_OmegaOmega];
    b = [b_Psi;
         b_Omega];
    M = sparse(M); % Enforce Sparse
    sol = M\b;
    Psi = Psi + sol(1:i_max*j_max);
    Omega = Omega + sol((i_max*j_max)+1:2*i_max*j_max);

    % Compute the new residual
    residual = norm(sol)/(2.0*i_max*j_max); % Normalized for # of equations being solved and DOF

    tTot = tTot + toc(tStart);

    % Plot solution
    if mod(iter,1) == 0
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
fprintf(fid, 'Total CPU-time: %s s\nAverage Time per Iteration: %s s\nPseudotime Step: %s', string(tTot), string(tAvg), string(Deltatau));
fclose(fid);

%% Functions
% ------------------------------------------------------------------------------

function k = pmap(i, j, i_max)
    k = i + (j-1)*i_max;
end