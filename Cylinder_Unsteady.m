%% Cylinder_Generalized
% ------------------------------------------------------------------------------
% This is a fully implicit solver for laminar NS Equations given in curvilinear
% coordinates as defined by generalized_o_grid_metrics for a cylinder with
% O-grid topology
% ------------------------------------------------------------------------------
clear all; close all; clc;
%
%% Computational Parameters
% ------------------------------------------------------------------------------
%
i_max = 181;
j_max = 181;
n_max = i_max*j_max;
%
% Inner and outer radii for cylinder and far-field
inner_rad = 0.5;
outer_rad = 50.0;
%
% Model 2D, Incompressible, Laminar, NS equations in Generalized Coordinates
%
% Flow conditions
Re = 40.0;
u_theta = 0.0;
u_infty = 1.0;
%
% The scheme should be stable for any timestep but in many cases if the
% initial guess/condition is far from the converged solution (think about
% this being similar to a Newton's method) the dtau value may need to be
% limited. A finite dtau effectively corresponds to underrelaxation in
% Newton's solution.
dtau = 1000.0;
% And double-timestepping
dt = 1.0;
nstart = 1;
%
% % Define parameters for Newton iteration
% residual = 1.0E5; % init residual
% epsilon = 1.0E-16; % drive residual down to this value before terminating
%
% File Info
mydir='C:\\Users\\Bucky\\Downloads\\cylinderSS_Results_Unsteady';
subfolder='Re'+string(Re)+'_alpha'+string(abs(u_theta))+'_ri'+string(inner_rad);
mkdir(fullfile(mydir,subfolder));
%
%% Grid Metrics
% ------------------------------------------------------------------------------
% This code sets up the metrics for 
% laminar NS Equations given in curvilinear
% coordinates as defined/given by Garmann (AIAA J. Vol. 51, No: 9)
%
% The current pointer system is for a cylinder with O-grid topology. 
% However, it is straightforward to modify the computation of metrics to 
% other topologies like C-grids, etc.
% ------------------------------------------------------------------------------
%
% Set dksi and deta
dksi = 1.0/(i_max-1);
deta = 1.0/(j_max-1);
%
% Arrays for pointers
kc   = zeros(i_max,j_max);
ke   = zeros(i_max,j_max);
kw   = zeros(i_max,j_max);
ks   = zeros(i_max,j_max);
kn   = zeros(i_max,j_max);
kne  = zeros(i_max,j_max);
knw  = zeros(i_max,j_max);
kse  = zeros(i_max,j_max);
ksw  = zeros(i_max,j_max);
%
% Initialize arrays for metrics
eta    = zeros(i_max*j_max,1);
ksi    = zeros(i_max*j_max,1);
x      = zeros(i_max*j_max,1);
y      = zeros(i_max*j_max,1);
Jac    = zeros(i_max*j_max,1);
alfa   = zeros(i_max*j_max,1);
beta   = zeros(i_max*j_max,1);
gama   = zeros(i_max*j_max,1);
P      = zeros(i_max*j_max,1);
Q      = zeros(i_max*j_max,1);
%
% Form the pointer system
for i = 1:i_max
    for j = 1:j_max
        kc(i,j) = i_max*(j-1) + i;
        ke(i,j) = kc(i,j) + 1;
        kw(i,j) = kc(i,j) - 1;
        kn(i,j) = kc(i,j) + i_max;
        ks(i,j) = kc(i,j) - i_max;
        %
        if i == 1
            kw(i,j) = kc(i,j) - 2 + i_max;
        elseif i == i_max
            ke(i,j) = kc(i,j) + 2 - i_max;
        end
        %
        kne(i,j) = ke(i,j) + i_max;
        knw(i,j) = kw(i,j) + i_max;
        kse(i,j) = ke(i,j) - i_max;
        ksw(i,j) = kw(i,j) - i_max;
        %
        % Define mesh in computational domain
        ksi(kc(i,j),1) = dksi*(i-1);
        eta(kc(i,j),1) = deta*(j-1);
    end
end
%
% Generate the computational grid using cosine clustering near the
% boundaries. This will allow better resolution of the boundary layer. It
% will also reduce the number of grid nodes to model the problem.
%
% Define the cylinder mesh in physical domain.     %
alen    = zeros(j_max,1);
alen    = 0.0;
alen(2) = 1.0;
alen(3) = 2.0;
alen(4) = 3.0;
for j = 5:j_max
    alen(j) = alen(j-1) + (j-3)^1;
end
alen2 = zeros(i_max,1);
alen2 = 0.0;
for i = 2:i_max
    alen2(i) = alen2(i-1) + min(i-1,i_max-i+1)^0.6;
end
alen  = alen/alen(j_max);
alen2 = alen2/alen2(i_max);
%
for i = 1:i_max
    for j = 1:j_max
        theta = 2*pi*(1 - alen2(i));
        rad   = inner_rad + (outer_rad - inner_rad)*alen(j);
        x(kc(i,j),1) = rad*cos(theta);
        y(kc(i,j),1) = rad*sin(theta);
    end
end
%
% Map the 1-D grid arrays to 2-D for plotting
for i = 1:i_max
    for j = 1:j_max
        xg(i,j) = x(kc(i,j),1);
        yg(i,j) = y(kc(i,j),1);
    end
end
%
%Plot the computational mesh
figure(1);
mesh(xg,yg,0*xg);
axis([-2 16 -4 4]);
axis equal;
axis([-2 16 -4 4]);
view(0,90);
drawnow;
% pause;
%
% First derivatives
dxdksi = zeros(i_max*j_max,1);
dxdeta = zeros(i_max*j_max,1);
dydksi = zeros(i_max*j_max,1);
dydeta = zeros(i_max*j_max,1);
%
% Second derivatives
ddxdksidksi = zeros(i_max*j_max,1);
ddxdetadeta = zeros(i_max*j_max,1);
ddxdksideta = zeros(i_max*j_max,1);
ddydksidksi = zeros(i_max*j_max,1);
ddydetadeta = zeros(i_max*j_max,1);
ddydksideta = zeros(i_max*j_max,1);
%
%
% Compute dxdksi, dydksi, dxdeta, dydeta, etc.
% Note that (x,y) and (xi,eta) are related to each other
% Inner points first - treat i = 1 and i = i_max as inner since they
% correspond to the branch cut
for i = 1:i_max
    for j = 2:j_max-1
        %
        % First derivatives
        dxdksi(kc(i,j),1) = (x(ke(i,j),1) - x(kw(i,j),1))/2/dksi;
        dydksi(kc(i,j),1) = (y(ke(i,j),1) - y(kw(i,j),1))/2/dksi;
        dxdeta(kc(i,j),1) = (x(kn(i,j),1) - x(ks(i,j),1))/2/deta;
        dydeta(kc(i,j),1) = (y(kn(i,j),1) - y(ks(i,j),1))/2/deta;
        %
        % Second derivatives
        ddxdksidksi(kc(i,j),1) = (x(ke(i,j),1) - 2*x(kc(i,j),1) + x(kw(i,j),1))/dksi^2;
        ddydksidksi(kc(i,j),1) = (y(ke(i,j),1) - 2*y(kc(i,j),1) + y(kw(i,j),1))/dksi^2;
        ddxdetadeta(kc(i,j),1) = (x(kn(i,j),1) - 2*x(kc(i,j),1) + x(ks(i,j),1))/deta^2;
        ddydetadeta(kc(i,j),1) = (y(kn(i,j),1) - 2*y(kc(i,j),1) + y(ks(i,j),1))/deta^2;
    end
end
%
% Bottom Boundary (xi-derivatives)
for i = 1:i_max
    for j = 1:1
        dxdksi(kc(i,j),1) = (x(ke(i,j),1) - x(kw(i,j),1))/2/dksi;
        dydksi(kc(i,j),1) = (y(ke(i,j),1) - y(kw(i,j),1))/2/dksi;
        %
        ddxdksidksi(kc(i,j),1) = (x(ke(i,j),1) - 2*x(kc(i,j),1) + x(kw(i,j),1))/dksi^2;
        ddydksidksi(kc(i,j),1) = (y(ke(i,j),1) - 2*y(kc(i,j),1) + y(kw(i,j),1))/dksi^2;
    end
end
%
% Bottom Boundary (eta-derivatives)
for i = 1:i_max
    for j = 1:1
        dxdeta(kc(i,j),1) = (-3*x(kc(i,j),1) + 4*x(kn(i,j),1) - x(kn(i,j+1),1))/2/deta;
        dydeta(kc(i,j),1) = (-3*y(kc(i,j),1) + 4*y(kn(i,j),1) - y(kn(i,j+1),1))/2/deta;
        %
        ddxdetadeta(kc(i,j),1) = (2*x(kc(i,j),1) - 5*x(kn(i,j),1) + 4*x(kn(i,j+1),1) - x(kn(i,j+2),1))/deta^2;
        ddydetadeta(kc(i,j),1) = (2*y(kc(i,j),1) - 5*y(kn(i,j),1) + 4*y(kn(i,j+1),1) - y(kn(i,j+2),1))/deta^2;
    end
end
%
% Top Boundary (xi-derivatives)
for i = 1:i_max
    for j = j_max:j_max
        dxdksi(kc(i,j),1) = (x(ke(i,j),1) - x(kw(i,j),1))/2/dksi;
        dydksi(kc(i,j),1) = (y(ke(i,j),1) - y(kw(i,j),1))/2/dksi;
        %
        ddxdksidksi(kc(i,j),1) = (x(ke(i,j),1) - 2*x(kc(i,j),1) + x(kw(i,j),1))/dksi^2;
        ddydksidksi(kc(i,j),1) = (y(ke(i,j),1) - 2*y(kc(i,j),1) + y(kw(i,j),1))/dksi^2;
    end
end
%
% Top Boundary (eta-derivatives)
for i = 1:i_max
    for j = j_max:j_max
        dxdeta(kc(i,j),1) = ( 3*x(kc(i,j),1) - 4*x(ks(i,j),1) + x(ks(i,j-1),1))/2/deta;
        dydeta(kc(i,j),1) = ( 3*y(kc(i,j),1) - 4*y(ks(i,j),1) + y(ks(i,j-1),1))/2/deta;
        %
        ddxdetadeta(kc(i,j),1) = (2*x(kc(i,j),1) - 5*x(ks(i,j),1) + 4*x(ks(i,j-1),1) - x(ks(i,j-2),1))/deta^2;
        ddydetadeta(kc(i,j),1) = (2*y(kc(i,j),1) - 5*y(ks(i,j),1) + 4*y(ks(i,j-1),1) - y(ks(i,j-2),1))/deta^2;
    end
end
%
% Cross-derivatives everywhere
for i = 1:i_max
    for j = 1:j_max
        ddxdksideta(kc(i,j),1) = (dxdeta(ke(i,j),1) - dxdeta(kw(i,j),1))/2/dksi;
        ddydksideta(kc(i,j),1) = (dydeta(ke(i,j),1) - dydeta(kw(i,j),1))/2/dksi;
    end
end
%
% Compute the Jacobian of transformation using physical coordinates (x,y)
for i = 1:i_max
    for j = 1:j_max
        Jac(kc(i,j),1)  = 1.0/(dxdksi(kc(i,j),1)*dydeta(kc(i,j),1) - dxdeta(kc(i,j),1)*dydksi(kc(i,j),1));
    end
end
%
% Now do the transformation metrics for the computational coor. (xi,eta)
% dksidx, dksidy, detadx, detady as well as the coefficients, alfa, beta,
% gama, P, and Q as defined by Garmann
for i = 1:i_max
    for j = 1:j_max
        k = kc(i,j);
        dksidx(k,1) = Jac(k,1)*dydeta(k,1);
        dksidy(k,1) =-Jac(k,1)*dxdeta(k,1);
        detadx(k,1) =-Jac(k,1)*dydksi(k,1);
        detady(k,1) = Jac(k,1)*dxdksi(k,1);
        %
        alfa(k,1) = dksidx(k,1)^2 + dksidy(k,1)^2;
        beta(k,1) = detadx(k,1)^2 + detady(k,1)^2;
        gama(k,1) = dksidx(k,1)*detadx(k) + dksidy(k)*detady(k);
        %
        P(k,1)    =-( 1*alfa(k,1)*(ddxdksidksi(k,1)*dksidx(k,1) + ddydksidksi(k,1)*dksidy(k,1))...
                     +2*gama(k,1)*(ddxdksideta(k,1)*dksidx(k,1) + ddydksideta(k,1)*dksidy(k,1))...
                     +1*beta(k,1)*(ddxdetadeta(k,1)*dksidx(k,1) + ddydetadeta(k,1)*dksidy(k,1)));
        %       
        Q(k,1)    =-( 1*alfa(k,1)*(ddxdksidksi(k,1)*detadx(k,1) + ddydksidksi(k,1)*detady(k,1))...
                     +2*gama(k,1)*(ddxdksideta(k,1)*detadx(k,1) + ddydksideta(k,1)*detady(k,1))...
                     +1*beta(k,1)*(ddxdetadeta(k,1)*detadx(k,1) + ddydetadeta(k,1)*detady(k,1)));
    end
end
%
%% CFD Script
% ------------------------------------------------------------------------------
% Discretize and convert to a linear system A*Δx = b
%
% Init
A_PsiPsi = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Allocate Streamfxn Coeff Matrix Sparsely
A_PsiOmega = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Allocate Coeff Matrix for Streamfxn Dependence on Vorticity
A_OmegaOmega = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Allocate Vorticity Coeff Matrix Sparsely
A_OmegaPsi = spalloc(i_max*j_max, i_max*j_max, 5*i_max*j_max); % Allocate Coeff Matrix for Vorticity Dependence on Streamfxn
b_Psi = zeros(i_max*j_max,1); % init RHS
b_Omega = zeros(i_max*j_max,1);
%
% Define Solution Variables (1D because we use pointer mapping)
Psi = zeros(i_max*j_max,1);
Omega = zeros(i_max*j_max,1);
OmegaOld = zeros(i_max*j_max,1);
OmegaOldOld = zeros(i_max*j_max,1);
u = zeros(i_max*j_max,1);
v = zeros(i_max*j_max,1);
%
% Start Solver
% tic
nend = 10001;
for n = nstart:nend
    if n == 1
        time(n) = 0.0;
        dtdum = dt;
        dt = inf;
        dtaudum = dtau;
        dtau = inf;
    else
        dt = dtdum;
        time(n) = (n-1)*dt;
        dtau = inf;
    end
    if n >= 2 && n <= 3
        u_theta = 1.5;
    else
        u_theta = 0.0;
    end
    residual = 1.0E5; % reinit residual
    epsilon = 1.0E-16; % redefine epsilon -- drive residual down to this value before terminating
    iter = 0; % reinit iteration count
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
            % A_OmegaOmega(kc(i,j),kc(i,j)) = 1.0/dtau + (2.0*alfa(kc(i,j),1)/dksi^2 + 2.0*beta(kc(i,j),1)/deta^2)/Re;
            A_OmegaOmega(kc(i,j),kc(i,j)) = 3.0/(2.0*dt) + 1.0/dtau + (2.0*alfa(kc(i,j),1)/dksi^2 + 2.0*beta(kc(i,j),1)/deta^2)/Re;
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
    % Update
    OmegaOldOld = OmegaOld;
    OmegaOld = Omega;
    %
    % Begin iteration
    tTot = 0;
    iter = 0;
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
                                    Q(kc(i,j),1)*(Omega(kn(i,j),1)-Omega(ks(i,j),1))/(2.0*deta))/Re + ...
                                    -((3.0*Omega(kc(i,j),1)-4.0*OmegaOld(kc(i,j),1)+OmegaOldOld(kc(i,j),1))/(2.0*dt));
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
    %% Output Results
    % ------------------------------------------------------------------------------
    %
    % Stats on this Iteration
    fprintf(1, 'Time: %g, Total CPU-time: %s\n', time(n), string(tTot));
    %
    % Store an Arbitrary State Value in the Wake to Assess Shedding Frequency
    % I've chosen 10 degrees off of branch cut, 1/10 along branch cut
    myStateVar(n) = Omega(kc(floorDiv(i_max,36), floorDiv(j_max,10)));
    %
    % Plot Results
    figure(1);
    subplot(3, 1, 1);
    % Plot level curves for vorticity - arbitrary levels
    contour(xg,yg,reshape(Omega, i_max, j_max),[-6.0 -5.0 -4.0 -3.0 -2.0 -1.0 -0.5 -0.2 -0.1 -0.05 0.05 0.1 0.2 0.5 1.0 2.0 3.0 4.0 5.0 6.0],'LineWidth',2.0);
    axis([-2 16 -4 4]);
    axis equal;
    axis([-2 16 -4 4]);
    ylabel('y');
    xlabel('x');
    title('Vorticity Contour');
    subplot(3, 1, 2);
    % Plot level curve for streamfxn - levels informed by Ingham
    contour(xg,yg,reshape(Psi, i_max, j_max),[-0.5 -0.3 -0.25 -0.2 -0.15 -0.1 -0.08 -0.05 -0.04 -0.02 -0.01 0.01 0.02 0.04 0.05 0.08 0.1 0.15 0.2 0.25 0.3 0.5], 'LineWidth',2.0)
    axis([-2 16 -4 4]);
    axis equal;
    axis([-2 16 -4 4]);
    ylabel('y');
    xlabel('x');
    title('Streamline Pattern');
    shedding=subplot(3, 1, 3);
    % Shedding Frequency
    plot(time,myStateVar)
    % xlim([0,dt*nend])
    ylim([0.0,0.01])
    ylabel('Vorticity Value in Arbitrary Wake Position');
    xlabel('Time (s)');
    title('Vortex Shedding Pattern in Wake');
    drawnow;
end
%
% Save Final Shedding Plot
myAxes=findobj(shedding,'Type','Axes');
exportgraphics(myAxes,fullfile(mydir,subfolder,subfolder+'_sheddingPattern.jpg'));
%
% And Solution Matrices
save(fullfile(mydir,subfolder,subfolder+'_Psi.mat'), 'Psi')
save(fullfile(mydir,subfolder,subfolder+'_Omega.mat'), 'Omega')
%