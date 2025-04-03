%% CFD Script
% ------------------------------------------------------------------------------
% Redef Timestep
dtau = 5000.0;
%
% Define Coefficient Matrix
for i = 1:i_max
    for j = 2:j_max-1
        k = kc(i,j);
        % pointer mapping goes row-by-row to assemble Coeff. Matrix
        % Streamfxn
        A_PsiOmega(kc(i,j),kc(i,j)) = -1.0/Re;
        %
        A_PsiPsi(kc(i,j),kc(i,j)) = 1.0/dtau + (2.0*alfa(kc(i,j),1)/dksi^2 + 2.0*beta(kc(i,j),1)/deta^2)/Re;
        A_PsiPsi(kc(i,j),ke(i,j)) = (-alfa(kc(i,j),1)/dksi^2 - P(kc(i,j),1)/(2.0*dksi))/Re;
        A_PsiPsi(kc(i,j),kw(i,j)) = (-alfa(kc(i,j),1)/dksi^2 + P(kc(i,j),1)/(2.0*dksi))/Re;
        A_PsiPsi(kc(i,j),kn(i,j)) = (-beta(kc(i,j),1)/deta^2 - Q(kc(i,j),1)/(2.0*deta))/Re;
        A_PsiPsi(kc(i,j),ks(i,j)) = (-beta(kc(i,j),1)/deta^2 + Q(kc(i,j),1)/(2.0*deta))/Re;
        A_PsiPsi(kc(i,j),kne(i,j)) = (-gama(kc(i,j,1))/(2.0*deta*dksi))/Re;
        A_PsiPsi(kc(i,j),knw(i,j)) = (gama(kc(i,j,1))/(2.0*deta*dksi))/Re;
        A_PsiPsi(kc(i,j),kse(i,j)) = (gama(kc(i,j,1))/(2.0*deta*dksi))/Re;
        A_PsiPsi(kc(i,j),ksw(i,j)) = (-gama(kc(i,j,1))/(2.0*deta*dksi))/Re;
        % Transport
        A_OmegaOmega(kc(i,j),kc(i,j)) = 1.0/dtau + (2.0*alfa(kc(i,j),1)/dksi^2 + 2.0*beta(kc(i,j),1)/deta^2)/Re;
        % A_OmegaOmega(kc(i,j),ke(i,j)) = xxx; % Non-Constant Terms
        % A_OmegaOmega(kc(i,j),kw(i,j)) = xxx;
        % A_OmegaOmega(kc(i,j),kn(i,j)) = xxx;
        % A_OmegaOmega(kc(i,j),ks(i,j)) = xxx;
        A_OmegaOmega(kc(i,j),kne(i,j)) = (-gama(kc(i,j,1))/(2.0*deta*dksi))/Re;
        A_OmegaOmega(kc(i,j),knw(i,j)) = (gama(kc(i,j,1))/(2.0*deta*dksi))/Re;
        A_OmegaOmega(kc(i,j),kse(i,j)) = (gama(kc(i,j,1))/(2.0*deta*dksi))/Re;
        A_OmegaOmega(kc(i,j),ksw(i,j)) = (-gama(kc(i,j,1))/(2.0*deta*dksi))/Re;
    end
end
A_PsiPsi = sparse(A_PsiPsi); % Enforce Coeffs sparse
A_PsiOmega = sparse(A_PsiOmega); 
A_OmegaOmega = sparse(A_OmegaOmega); 
A_OmegaPsi = sparse(A_OmegaPsi);
%
% Begin iteration
tTot = 0;
while (residual > epsilon)
    % For CPU-time metrics
    tStart = tic;
    %
    for i = 1:i_max
        for j = 2:j_max-1
            % Define velocity (u,v) based on current Psi values
            u(kc(i,j),1) = dksidy(kc(i,j),1)*(Psi(ke(i,j),1)-Psi(kw(i,j),1))/(2.0*dksi) + detady(kc(i,j),1)*(Psi(kn(i,j),1)-Psi(ks(i,j),1))/(2.0*deta);
            v(kc(i,j),1) = -(dksidx(kc(i,j),1)*(Psi(ke(i,j),1)-Psi(kw(i,j),1))/(2.0*dksi) + detadx(kc(i,j),1)*(Psi(kn(i,j),1)-Psi(ks(i,j),1))/(2.0*deta));
            %
            % Update vorticity coefficient matrix
            A_OmegaOmega(kc(i,j),ke(i,j)) = (Jac(kc(i,j),1)*(Psi(kn(i,j),1)-Psi(ks(i,j),1))/(4.0*deta*dksi)) + (-alfa(kc(i,j),1)/dksi^2 - P(kc(i,j),1)/(2.0*dksi))/Re;
            A_OmegaOmega(kc(i,j),kw(i,j)) = (-Jac(kc(i,j),1)*(Psi(kn(i,j),1)-Psi(ks(i,j),1))/(4.0*deta*dksi)) + (-alfa(kc(i,j),1)/dksi^2 + P(kc(i,j),1)/(2.0*dksi))/Re;
            A_OmegaOmega(kc(i,j),kn(i,j)) = (-Jac(kc(i,j),1)*(Psi(ke(i,j),1)-Psi(kw(i,j),1))/(4.0*deta*dksi)) + (-beta(kc(i,j),1)/deta^2 - Q(kc(i,j),1)/(2.0*deta))/Re;
            A_OmegaOmega(kc(i,j),ks(i,j)) = (Jac(kc(i,j),1)*(Psi(ke(i,j),1)-Psi(kw(i,j),1))/(4.0*deta*dksi)) + (-beta(kc(i,j),1)/deta^2 + Q(kc(i,j),1)/(2.0*deta))/Re;
            A_OmegaPsi(kc(i,j),ke(i,j)) = -Jac(kc(i,j),1)*(Omega(kn(i,j),1)-Omega(ks(i,j),1))/(4*deta*dksi);
            A_OmegaPsi(kc(i,j),kw(i,j)) = Jac(kc(i,j),1)*(Omega(kn(i,j),1)-Omega(ks(i,j),1))/(4*deta*dksi);
            A_OmegaPsi(kc(i,j),kn(i,j)) = Jac(kc(i,j),1)*(Omega(ke(i,j),1)-Omega(kw(i,j),1))/(4*deta*dksi);
            A_OmegaPsi(kc(i,j),ks(i,j)) = -Jac(kc(i,j),1)*(Omega(ke(i,j),1)-Omega(kw(i,j),1))/(4*deta*dksi);
            %
            % And RHS
            b_Omega(kc(i,j),1) = -Jac(kc(i,j),1)*((Psi(kn(i,j),1)-Psi(ks(i,j),1))*(Omega(ke(i,j),1)-Omega(kw(i,j),1))/(4.0*deta*dksi) - ...
                                (Psi(ke(i,j),1)-Psi(kw(i,j),1))*(Omega(kn(i,j),1)-Omega(ks(i,j),1))/(4.0*deta*dksi)) + ...
                                ( alfa(kc(i,j),1)*(Omega(ke(i,j),1)-2*Omega(kc(i,j),1)+Omega(kw(i,j),1))/dksi^2 + ...
                                2*gama(kc(i,j),1)*(Omega(kne(i,j),1)-Omega(knw(i,j),1)-Omega(kse(i,j),1)+Omega(ksw(i,j),1))/(4.0*deta*dksi) + ...
                                beta(kc(i,j),1)*(Omega(kn(i,j),1)-2*Omega(kc(i,j),1)+Omega(ks(i,j),1))/deta^2 + ...
                                P(kc(i,j),1)*(Omega(ke(i,j),1)-Omega(kw(i,j),1))/(2.0*dksi) + ...
                                Q(kc(i,j),1)*(Omega(kn(i,j),1)-Omega(ks(i,j),1))/(2.0*deta))/Re;
            b_Psi(kc(i,j),1) = (Omega(kc(i,j),1) + ...
                                alfa(kc(i,j),1)*(Psi(ke(i,j),1)-2*Psi(kc(i,j),1)+Psi(kw(i,j),1))/dksi^2 + ...
                                2*gama(kc(i,j),1)*(Psi(kne(i,j),1)-Psi(knw(i,j),1)-Psi(kse(i,j),1)+Psi(ksw(i,j),1))/(4.0*deta*dksi) + ...
                                beta(kc(i,j),1)*(Psi(kn(i,j),1)-2*Psi(kc(i,j),1)+Psi(ks(i,j),1))/deta^2 + ...
                                P(kc(i,j),1)*(Psi(ke(i,j),1)-Psi(kw(i,j),1))/(2.0*dksi) + ...
                                Q(kc(i,j),1)*(Psi(kn(i,j),1)-Psi(ks(i,j),1))/(2.0*deta))/Re;
        end
    end
    %
    % Apply BCs
    % Cylinder Surface BC (Bottom of Computational Domain)
    for i = 1:i_max
        for j = 1:1
            % pointer mapping
            knn = kn(i,j) + i_max;
            theta = atan2(y(kc(i,j),1), x(kc(i,j),1));
            % r = sqrt(x(kc(i,j),1)^2 + y(kc(i,j),1)^2);
            % No-Slip
            A_PsiPsi(kc(i,j),kc(i,j)) = 1.0;
            b_Psi(kc(i,j),1) = -Psi(kc(i,j),1);
            % Rotational Velocity u_theta
            A_OmegaOmega(kc(i,j),kc(i,j)) = -1.0/Re;
            A_OmegaPsi(kc(i,j),kc(i,j)) = (7.0*beta(kc(i,j),1)/(2.0*deta^2))/Re;
            A_OmegaPsi(kc(i,j),kn(i,j)) = (-4.0*beta(kc(i,j),1)/(deta^2))/Re;
            A_OmegaPsi(kc(i,j),knn) = (beta(kc(i,j),1)/(2.0*deta^2))/Re;
            b_Omega(kc(i,j),1) = (Omega(kc(i,j),1) + ...
                                    beta(kc(i,j),1)*(-7.0*Psi(kc(i,j),1)+8.0*Psi(kn(i,j),1)-Psi(knn,1))/(2*deta^2) + ...
                                    (Q(kc(i,j),1)-3.0*beta(kc(i,j),1)/deta)*(-u_theta)/(detadx(kc(i,j),1)*cos(theta)+detady(kc(i,j),1)*sin(theta)))/Re;
            % fprintf(1,'%g\n',detadx(kc(i,j),1))
            % fprintf(1,'%g\n',sin(theta))
            % fprintf(1,'%g\n',detady(kc(i,j),1))
            % fprintf(1,'%g\n',cos(theta))
        end
    end
    % Far-Field BC (Top of Computational Domain)
    for i = 1:i_max
        for j = j_max:j_max
            % Crossflow Velocity u_infty
            A_PsiPsi(kc(i,j),kc(i,j)) = 1.0;
            b_Psi(kc(i,j),1) = u_infty*y(kc(i,j),1) - Psi(kc(i,j),1);
            % Vorticity to zero
            A_OmegaOmega(kc(i,j),kc(i,j)) = 1.0;
            b_Omega(kc(i,j),1) = -Omega(kc(i,j),1);
        end
    end
    %
    % Solve the linear systems for Psi and Omega
    M = [A_PsiPsi, A_PsiOmega;
         A_OmegaPsi, A_OmegaOmega];
    b = [b_Psi;
         b_Omega];
    M = sparse(M); % Enforce Sparse
    sol = M\b;
    Psi = Psi + sol(1:i_max*j_max);
    Omega = Omega + sol((i_max*j_max)+1:2*i_max*j_max);
    %
    % Compute the new residual
    residual = norm(sol)/(2.0*i_max*j_max); % Normalized for # of equations being solved and DOF
    %
    % Complete CPU-time metric cycle
    tTot = tTot + toc(tStart);
    %
    fprintf(1,'iter = %i, residual = %g\n',iter,log10(residual));
    iter = iter+1;
end
%
%% Output Results
% ------------------------------------------------------------------------------
%
% Plot Results
% uplot = reshape(u, i_max, j_max);
% vplot = reshape(v, i_max, j_max);
figure(1);
% Plot level curves for vorticity - arbitrary levels
contour(xg,yg,reshape(Omega, i_max, j_max),[-6.0 -5.0 -4.0 -3.0 -2.0 -1.0 -0.5 -0.2 -0.1 -0.05 0.05 0.1 0.2 0.5 1.0 2.0 3.0 4.0 5.0 6.0],'LineWidth',2.0);
axis([-2 16 -4 4]);
axis equal;
axis([-2 16 -4 4]);
ylabel('y');
xlabel('x');
title('Vorticity Contour');
figure(2);
% Plot level curve for streamfxn - levels informed by Ingham
contour(xg,yg,reshape(Psi, i_max, j_max),[-0.5 -0.3 -0.25 -0.2 -0.15 -0.1 -0.08 -0.05 -0.04 -0.02 -0.01 0.01 0.02 0.04 0.05 0.08 0.1 0.15 0.2 0.25 0.3 0.5], 'LineWidth',2.0)
axis([-2 16 -4 4]);
axis equal;
axis([-2 16 -4 4]);
ylabel('y');
xlabel('x');
title('Streamline Pattern');
figure(3);
% Plot surface vorticity against theta
myTheta = atan2(y(kc(1:i_max,1),1), x(kc(1:i_max,1),1));
myTheta(myTheta<0)=myTheta(myTheta<0)+(2*pi);
myTheta=myTheta*180/pi;
plot(myTheta,-Omega(kc(1:i_max,1),1))
xlim([0,360])
% plot(-Omega(kc(i_max:-1:1),1))
ylabel('Vorticity');
xlabel('Angle');
title('The Variation of Surface Vorticity');
drawnow;
%
% Save Final Plots
saveas(figure(1),fullfile(mydir,subfolder,subfolder+'_vorticity.jpg'));
saveas(figure(2),fullfile(mydir,subfolder,subfolder+'_streamfxn.jpg'));
saveas(figure(3),fullfile(mydir,subfolder,subfolder+'_surfaceVorticity.jpg'));
%
% And Solution Matrices
save(fullfile(mydir,subfolder,subfolder+'_Psi.mat'), 'Psi')
save(fullfile(mydir,subfolder,subfolder+'_Omega.mat'), 'Omega')
%
% And Timing Info
tAvg = tTot/iter;
fid = fopen(fullfile(mydir,subfolder,'time.txt'),'wt');
fprintf(fid, 'Total CPU-time: %s s\nAverage Time per Iteration: %s s\nPseudotime Step: %s', string(tTot), string(tAvg), string(dtau));
fclose(fid);