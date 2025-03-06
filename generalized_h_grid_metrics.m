% This code sets up the metrics for 
% laminar NS Equations given in curvilinear
% coordinates as defined/given by Garmann (AIAA J. Vol. 51, No: 9)
%
% The current pointer system is for a square cavity. However, it is
% straightforward to modify the computation of metrics to model flow around
% objects such as airfoils, cylinders, etc.
clear all;
close all;
%
%imax = input(' The number of points in  xi?\n');
%jmax = input(' The number of points in eta?\n');
%
imax = 101;
jmax = 101;
nmax = imax*jmax;
%
% Metrics for NS equations in Generalized Coordinates
%
% dksi and deta
%dksi = 1.0;
%deta = 1.0;
dksi = 1.0/(imax-1);
deta = 1.0/(jmax-1);
%
dx = 1.0/(imax-1);
dy = 1.0/(jmax-1);
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
Jac2   = zeros(imax*jmax,1);
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
% Define the cavity mesh in physical domain. Note: the mesh is clustered
% near the boundaries using a cosine distribution. alen is an array with
% normalized values of arclength, i.e, alen E [0 1] that provides clustering.
alen = zeros(imax,1);
alen = 0.0;
%
% Cosine distribution
% for i = 2:imax
%    alen(i) = 0.5 - 0.5*cos(pi*(i-1)/(imax-1)); % this is for cosine clustering
% end
%
% Equal spacing
for i = 2:imax
     alen(i) = alen(i-1) + 1; % This is for equally spaced mesh.
end
%
% Algebraic
for i = 2:(imax+1)/2
alen(i) = alen(i-1) + (i-1)^1; % This is for geometric clustering.
end
for i = (imax+1)/2+1:imax
alen(i) = alen(i-1) + (imax-i+1)^1; % This is for geometrix clustering.
end
%
alen = alen/alen(imax);
%
for i = 1:imax
    for j = 1:jmax
        x(kc(i,j),1) = alen(i);
        y(kc(i,j),1) = alen(j);
    end
end
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
% Inner points first
for i = 2:imax-1
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
        ddxdksideta(kc(i,j),1) = (x(kne(i,j),1) - x(kse(i,j),1) - x(knw(i,j),1) + x(ksw(i,j),1))/4/dksi/deta;
        ddydksideta(kc(i,j),1) = (y(kne(i,j),1) - y(kse(i,j),1) - y(knw(i,j),1) + y(ksw(i,j),1))/4/dksi/deta;
    end
end
%
% Left Boundary (xi-derivatives)
for i = 1:1
    for j = 1:jmax
        dxdksi(kc(i,j),1) = (-3*x(kc(i,j),1) + 4*x(ke(i,j),1) - x(ke(i+1,j),1))/2/dksi;
        dydksi(kc(i,j),1) = (-3*y(kc(i,j),1) + 4*y(ke(i,j),1) - y(ke(i+1,j),1))/2/dksi;
        %
        ddxdksidksi(kc(i,j),1) = (2*x(kc(i,j),1) - 5*x(ke(i,j),1) + 4*x(ke(i+1,j),1) - x(ke(i+2,j),1))/dksi^2;
        ddydksidksi(kc(i,j),1) = (2*y(kc(i,j),1) - 5*y(ke(i,j),1) + 4*y(ke(i+1,j),1) - y(ke(i+2,j),1))/dksi^2;
    end
end
%
% Left Boundary (eta-derivatives)
for i = 1:1
    for j = 2:jmax-1
        dxdeta(kc(i,j),1) = (x(kn(i,j),1) - x(ks(i,j),1))/2/deta;
        dydeta(kc(i,j),1) = (y(kn(i,j),1) - y(ks(i,j),1))/2/deta;
        %
        ddxdetadeta(kc(i,j),1) = (x(kn(i,j),1) - 2*x(kc(i,j),1) + x(ks(i,j),1))/deta^2;
        ddydetadeta(kc(i,j),1) = (y(kn(i,j),1) - 2*y(kc(i,j),1) + y(ks(i,j),1))/deta^2;
    end
end
%
% Right Boundary (xi-derivatives)
for i = imax:imax
    for j = 1:jmax
        dxdksi(kc(i,j),1) = ( 3*x(kc(i,j),1) - 4*x(kw(i,j),1) + x(kw(i-1,j),1))/2/dksi;
        dydksi(kc(i,j),1) = ( 3*y(kc(i,j),1) - 4*y(kw(i,j),1) + y(kw(i-1,j),1))/2/dksi;
        %
        ddxdksidksi(kc(i,j),1) = (2*x(kc(i,j),1) - 5*x(kw(i,j),1) + 4*x(kw(i-1,j),1) - x(kw(i-2,j),1))/dksi^2;
        ddydksidksi(kc(i,j),1) = (2*y(kc(i,j),1) - 5*y(kw(i,j),1) + 4*y(kw(i-1,j),1) - y(kw(i-2,j),1))/dksi^2;
    end
end
%
% Right Boundary (eta-derivatives)
for i = imax:imax
    for j = 2:jmax-1
        dxdeta(kc(i,j),1) = (x(kn(i,j),1) - x(ks(i,j),1))/2/deta;
        dydeta(kc(i,j),1) = (y(kn(i,j),1) - y(ks(i,j),1))/2/deta;
        %
        ddxdetadeta(kc(i,j),1) = (x(kn(i,j),1) - 2*x(kc(i,j),1) + x(ks(i,j),1))/deta^2;
        ddydetadeta(kc(i,j),1) = (y(kn(i,j),1) - 2*y(kc(i,j),1) + y(ks(i,j),1))/deta^2;
    end
end
%
% Bottom Boundary (xi-derivatives)
for i = 2:imax-1
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
for i = 2:imax-1
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
% Cross-derivatives on the top and bottom excluding i = 1 and i = imax; use eta derivatives
for i = 2:imax-1
    for j = 1:jmax-1:jmax
        ddxdksideta(kc(i,j),1) = (dxdeta(ke(i,j),1) - dxdeta(kw(i,j),1))/2/dksi;
        ddydksideta(kc(i,j),1) = (dydeta(ke(i,j),1) - dydeta(kw(i,j),1))/2/dksi;
    end
end
%
% Cross-derivatives at i = 1 and i = imax excluding the top and bottom; use xi derivatives
for i = 1:imax-1:imax
    for j = 2:jmax-1
        ddxdksideta(kc(i,j),1) = (dxdksi(kn(i,j),1) - dxdksi(ks(i,j),1))/2/deta;
        ddydksideta(kc(i,j),1) = (dydksi(kn(i,j),1) - dydksi(ks(i,j),1))/2/deta;
    end
end
%
% Cross-derivatives at 4 corners
i = 1;
j = 1;
ddxdksideta(kc(i,j),1) = (-3*dxdksi(kc(i,j),1) + 4*dxdksi(kn(i,j),1) - dxdksi(kn(i,j+1),1))/2/deta;
ddydksideta(kc(i,j),1) = (-3*dydksi(kc(i,j),1) + 4*dydksi(kn(i,j),1) - dydksi(kn(i,j+1),1))/2/deta;
i = 1;
j = jmax;
ddxdksideta(kc(i,j),1) = ( 3*dxdksi(kc(i,j),1) - 4*dxdksi(ks(i,j),1) + dxdksi(ks(i,j-1),1))/2/deta;
ddydksideta(kc(i,j),1) = ( 3*dydksi(kc(i,j),1) - 4*dydksi(ks(i,j),1) + dydksi(ks(i,j-1),1))/2/deta;
i = imax;
j = 1;
ddxdksideta(kc(i,j),1) = (-3*dxdksi(kc(i,j),1) + 4*dxdksi(kn(i,j),1) - dxdksi(kn(i,j+1),1))/2/deta;
ddydksideta(kc(i,j),1) = (-3*dydksi(kc(i,j),1) + 4*dydksi(kn(i,j),1) - dydksi(kn(i,j+1),1))/2/deta;
i = imax;
j = jmax;
ddxdksideta(kc(i,j),1) = ( 3*dxdksi(kc(i,j),1) - 4*dxdksi(ks(i,j),1) + dxdksi(ks(i,j-1),1))/2/deta;
ddydksideta(kc(i,j),1) = ( 3*dydksi(kc(i,j),1) - 4*dydksi(ks(i,j),1) + dydksi(ks(i,j-1),1))/2/deta;
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
% This is where the CFD solver starts
xg = reshape(x,imax,jmax);
yg = reshape(y,imax,jmax);

mesh(xg,yg,0*xg)
contour(xg,yg,reshape(Jac,imax,jmax))

