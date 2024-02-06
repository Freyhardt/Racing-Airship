
%% Author: Jan S. Freyhardt \n,
%% Bachelor Thesis \n,
%% Eight Funktions
%% email: jan.freyhardt@gmail.com\n

clc
close all
clear

t = linspace(0, 2*pi, 1000);  % Generates 1000 equidistant points in the interval [0, 2*pi]

%Funktions
x = sin(t);
y = sin(t) .* cos(t);

% Greating a Plot
plot(x, y);
axis equal; % Ensures that the axes are drawn in a 1:1 ratio
title('Eight Slope Equation');
xlabel('x');
ylabel('y');
