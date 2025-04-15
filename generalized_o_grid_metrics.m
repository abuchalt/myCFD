% This code sets up the metrics for 
% laminar NS Equations given in curvilinear
% coordinates as defined/given by Garmann (AIAA J. Vol. 51, No: 9)
%
% The current pointer system is for a cylinder with O-grid topology. 
% However, it is straightforward to modify the computation of metrics to 
% other topologies like C-grids, etc. 
clear all;
close all;
%
imax = 181;
jmax = 181;
nmax = imax*jmax;
%
% Inner and outer radii for cylinder and far-field
inner_rad = 0.50;
outer_rad = 50.0;
%
% Model 2D, Incompressible, Laminar, NS equations in Generalized Coordinates
%
% Flow conditions
Re   = 100.0;
%
% The scheme should be stable for any timestep but in many cases if the
% initial guess/condition is far from the converged solution (think about
% this being similar to a Newton's method) the dtau value may need to be
% limited. A finite dtau effectively corresponds to underrelaxation in
% Newton's solution.
dtau = inf;
%
% Set dksi and deta
dksi = 1.0/(imax-1);
deta = 1.0/(jmax-1);
%
% Arrays for pointers
kc   = zeros(imax,jmax);
ke   = zeros(imax,jmax);
kw   = zeros(imax,jmax);
ks   = zeros(imax,jmax);
kn   = zeros(imax,jmax);
kne  = zeros(imax,jmax);
knw  = zeros(imax,jmax);
kse  = zeros(imax,jmax);
ksw  = zeros(imax,jmax);
%
% Initialize arrays for metrics
eta    = zeros(imax*jmax,1);
ksi    = zeros(imax*jmax,1);
x      = zeros(imax*jmax,1);
y      = zeros(imax*jmax,1);
Jac    = zeros(imax*jmax,1);
alfa   = zeros(imax*jmax,1);
beta   = zeros(imax*jmax,1);
gama   = zeros(imax*jmax,1);
P      = zeros(imax*jmax,1);
Q      = zeros(imax*jmax,1);
%
% Form the pointer system
for i = 1:imax
    for j = 1:jmax
        kc(i,j) = imax*(j-1) + i;
        ke(i,j) = kc(i,j) + 1;
        kw(i,j) = kc(i,j) - 1;
        kn(i,j) = kc(i,j) + imax;
        ks(i,j) = kc(i,j) - imax;
        %
        if i == 1
            kw(i,j) = kc(i,j) - 2 + imax;
        elseif i == imax
            ke(i,j) = kc(i,j) + 2 - imax;
        end
        %
        kne(i,j) = ke(i,j) + imax;
        knw(i,j) = kw(i,j) + imax;
        kse(i,j) = ke(i,j) - imax;
        ksw(i,j) = kw(i,j) - imax;
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
alen    = zeros(jmax,1);
% alen    = 0.0;
alen(2) = 1.0;
alen(3) = 2.0;
alen(4) = 3.0;
for j = 5:jmax
    alen(j) = alen(j-1) + (j-3)^1;
end
alen2 = zeros(imax,1);
% alen2 = 0.0;
for i = 2:imax
    alen2(i) = alen2(i-1) + min(i-1,imax-i+1)^0.6;
end
alen  = alen/alen(jmax);
alen2 = alen2/alen2(imax);
%
for i = 1:imax
    for j = 1:jmax
        theta = 2*pi*(1 - alen2(i));
        rad   = inner_rad + (outer_rad - inner_rad)*alen(j);
        x(kc(i,j),1) = rad*cos(theta);
        y(kc(i,j),1) = rad*sin(theta);
    end
end
%
% Map the 1-D grid arrays to 2-D for plotting
for i = 1:imax
    for j = 1:jmax
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
pause;
%
% First derivatives
dxdksi = zeros(imax*jmax,1);
dxdeta = zeros(imax*jmax,1);
dydksi = zeros(imax*jmax,1);
dydeta = zeros(imax*jmax,1);
%
% Second derivatives
ddxdksidksi = zeros(imax*jmax,1);
ddxdetadeta = zeros(imax*jmax,1);
ddxdksideta = zeros(imax*jmax,1);
ddydksidksi = zeros(imax*jmax,1);
ddydetadeta = zeros(imax*jmax,1);
ddydksideta = zeros(imax*jmax,1);
%
%
% Compute dxdksi, dydksi, dxdeta, dydeta, etc.
% Note that (x,y) and (xi,eta) are related to each other
% Inner points first - treat i = 1 and i = imax as inner since they
% correspond to the branch cut
for i = 1:imax
    for j = 2:jmax-1
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
for i = 1:imax
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
for i = 1:imax
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
for i = 1:imax
    for j = jmax:jmax
        dxdksi(kc(i,j),1) = (x(ke(i,j),1) - x(kw(i,j),1))/2/dksi;
        dydksi(kc(i,j),1) = (y(ke(i,j),1) - y(kw(i,j),1))/2/dksi;
        %
        ddxdksidksi(kc(i,j),1) = (x(ke(i,j),1) - 2*x(kc(i,j),1) + x(kw(i,j),1))/dksi^2;
        ddydksidksi(kc(i,j),1) = (y(ke(i,j),1) - 2*y(kc(i,j),1) + y(kw(i,j),1))/dksi^2;
    end
end
%
% Top Boundary (eta-derivatives)
for i = 1:imax
    for j = jmax:jmax
        dxdeta(kc(i,j),1) = ( 3*x(kc(i,j),1) - 4*x(ks(i,j),1) + x(ks(i,j-1),1))/2/deta;
        dydeta(kc(i,j),1) = ( 3*y(kc(i,j),1) - 4*y(ks(i,j),1) + y(ks(i,j-1),1))/2/deta;
        %
        ddxdetadeta(kc(i,j),1) = (2*x(kc(i,j),1) - 5*x(ks(i,j),1) + 4*x(ks(i,j-1),1) - x(ks(i,j-2),1))/deta^2;
        ddydetadeta(kc(i,j),1) = (2*y(kc(i,j),1) - 5*y(ks(i,j),1) + 4*y(ks(i,j-1),1) - y(ks(i,j-2),1))/deta^2;
    end
end
%
% Cross-derivatives everywhere
for i = 1:imax
    for j = 1:jmax
        ddxdksideta(kc(i,j),1) = (dxdeta(ke(i,j),1) - dxdeta(kw(i,j),1))/2/dksi;
        ddydksideta(kc(i,j),1) = (dydeta(ke(i,j),1) - dydeta(kw(i,j),1))/2/dksi;
    end
end
%
% Compute the Jacobian of transformation using physical coordinates (x,y)
for i = 1:imax
    for j = 1:jmax
        Jac(kc(i,j),1)  = 1.0/(dxdksi(kc(i,j),1)*dydeta(kc(i,j),1) - dxdeta(kc(i,j),1)*dydksi(kc(i,j),1));
    end
end
%
% Now do the transformation metrics for the computational coor. (xi,eta)
% dksidx, dksidy, detadx, detady as well as the coefficients, alfa, beta,
% gama, P, and Q as defined by Garmann
for i = 1:imax
    for j = 1:jmax
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
% Now you're ready for the CFD part

