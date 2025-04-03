% y' = -y^2; y(0) = 1
% analytical: y = 1/(t+1)
clear all; close all; clc
%
% Define RK4 inputs:
dt = 0.01;
y0 = 1.0;
%
for t = dt:dt:2.0
    %
    dy1 = dt*(-(y0)^2);
    dy2 = dt*(-(y0+dy1/2)^2);
    dy3 = dt*(-(y0+dy2/2)^2);
    dy4 = dt*(-(y0+dy3/1)^2);
    %+3d1
    y = y0 + 1/6*(dy1+2*dy2+2*dy3+dy4);
    %
    fprintf(1, 'time=%5.2f, y_exact = %12.8f, y = %12.8f, error=%e\n', t, 1/(1+t), y, abs(y-1/(1+t)));
    %
    y0 = y;
end
%
pause;
%
alpha = [1/4 1/3 1/2 1];
n_rk = 4;
for t = dt:dt:2.0
    y = y0;
    for n = 1:n_rk
        y = y + alpha(n)*dt*(-y^2);
    end
    %
    fprintf(1, 'time=%5.2f, y_exact = %12.8f, y = %12.8f, error=%e\n', t, 1/(1+t), y, abs(y-1/(1+t)));
    %
    y0 = y;
end
