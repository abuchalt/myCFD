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
% Flow conditions
Re = 100.0;
u_infty = 1.0;
D = inner_rad*2;
%
% And double-timestepping
dt = 0.05;
%
% File Info
mydir='C:\\Users\\Bucky\\Downloads\\cylinderSS_Results_Unsteady';
subfolder='Re'+string(Re)+'_unsteady_long';
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
% alen    = 0.0;
alen(2) = 1.0;
alen(3) = 2.0;
alen(4) = 3.0;
for j = 5:j_max
    alen(j) = alen(j-1) + (j-3)^1;
end
alen2 = zeros(i_max,1);
% alen2 = 0.0;
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
% Grab Matrix Files
load(fullfile(mydir,subfolder,subfolder+'_Psi.mat'), 'Psi')
load(fullfile(mydir,subfolder,subfolder+'_Omega.mat'), 'Omega')
%
load(fullfile(mydir,subfolder,subfolder+'_myStateVar.mat'), 'myStateVar')
load(fullfile(mydir,subfolder,subfolder+'_time.mat'), 'time')
%
% Plotting
figure(1);
% Plot level curves for vorticity - arbitrary levels
[M,c] = contour(xg,yg,reshape(Omega, i_max, j_max),[-6.0 -5.0 -4.0 -3.0 -2.0 -1.0 -0.5 -0.2 -0.1 -0.05 0.05 0.1 0.2 0.5 1.0 2.0 3.0 4.0 5.0 6.0],'LineWidth',2.0);
c.LineWidth = 0.5;
axis([-2 16 -4 4]);
axis equal;
axis([-2 16 -4 4]);
ylabel('y');
xlabel('x');
title('Vorticity Contour');
figure(2);
% Plot level curve for streamfxn - levels informed by Ingham
[M,c] = contour(xg,yg,reshape(Psi, i_max, j_max),[-0.5 -0.3 -0.25 -0.2 -0.15 -0.1 -0.08 -0.05 -0.04 -0.02 -0.01 0.01 0.02 0.04 0.05 0.08 0.1 0.15 0.2 0.25 0.3 0.5], 'LineWidth',2.0);
c.LineWidth = 0.5;
axis([-2 16 -4 4]);
axis equal;
axis([-2 16 -4 4]);
ylabel('y');
xlabel('x');
title('Streamline Pattern');
figure(3);
% Shedding Frequency
plot(time,myStateVar)
% xlim([0,dt*nend])
ylim([-0.6,0.8]);
pbaspect([2.5 1 1]);
ylabel('Vorticity Value in Arbitrary Wake Position');
xlabel('Time (unitless)');
title('Vortex Shedding Pattern in Wake');
figure(4);
freqplot = fft(myStateVar);
myN = length(freqplot);
fs = 1/dt;
f = (0:myN-1)*fs/myN;
fshift = (-myN/2:myN/2-1)*(fs/myN);
yshift = fftshift(freqplot);
% plot(f,abs(freqplot));
plot(fshift,abs(yshift))
xlabel('Frequency (inverse unitless time)');
ylabel('Magnitude');
title('Strouhal Number Determination');
drawnow;
[pks, locs] = findpeaks(abs(yshift), fshift); % Finds peaks and their locations
% Get the highest peak
[max_peak, max_index] = max(pks);
max_frequency = locs(max_index);
fprintf(1, 'Max Frequency: %g\n', max_frequency);
St = abs(max_frequency)*D/u_infty;
fprintf(1, 'Strouhal Number: %g\n', St);
%
% Save Final Plots
saveas(figure(1),fullfile(mydir,subfolder,subfolder+'_vorticity.jpg'));
saveas(figure(2),fullfile(mydir,subfolder,subfolder+'_streamfxn.jpg'));
saveas(figure(3),fullfile(mydir,subfolder,subfolder+'_sheddingPattern.jpg'));
saveas(figure(4),fullfile(mydir,subfolder,subfolder+'_strouhalFreq.jpg'));