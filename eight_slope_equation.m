
%% Author: Jan S. Freyhardt \n,
%% Bachelor Thesis \n,
%% Eight Slope Funktions with Pylons
%% email: jan.freyhardt@gmail.com\n

clc
clear
close all

% Parameters
a = 15;  % Maximum extension on the pylons
R_min = 9; % Minimum possible radius

% Parametric equation of the "8" with constant curvature radius
t = linspace(-pi/2, 3*pi/2, 1000); % Range corresponding to the modified equations
x = R_min*2 * sqrt(2) * cos(t) ./ (sin(t).^2 + 1);
y = R_min*2 * sqrt(2) * cos(t) .* sin(t) ./ (sin(t).^2 + 1);

pylon1 = [-a/2, 0];
pylon2 = [a/2, 0];

% Plot of the "8" with constant radius of curvature
figure;
plot(x, y, 'b', 'LineWidth', 2);
hold on;
plot(pylon1(1), pylon1(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(pylon2(1), pylon2(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Racing Line');
xlabel('x');
ylabel('y');
legend('Flight Line','Pylon 1','Pylon 2');