%% Author: Jan S. Freyhardt \n,
%% Bachelor Thesis \n,
%% Plane Optimisation
%% email: jan.freyhardt@gmail.com\n

clc
clear
close all

% Erstelle ein Gitter von Werten für S_p und S_r
S_p_values = linspace(0, 0.04574, 100); % Anpassen der Schritte und Grenzen nach Bedarf
S_r_values = linspace(0, 0.04574, 100);
[S_p_grid, S_r_grid] = meshgrid(S_p_values, S_r_values);

% Berechne den Funktionswert (dt) für jedes Paar von S_p und S_r
dt_values = zeros(size(S_p_grid));
for i = 1:numel(S_p_grid)
    dt_values(i) = compute_dt([S_p_grid(i), S_r_grid(i)]);
end

% Erstelle eine 3D-Oberflächengrafik
figure;
surf(S_p_grid, S_r_grid, dt_values);
xlabel('S_p');
ylabel('S_r');
zlabel('dt');
title('3D-Surface Graphic of the Function dt');

% Führe die Optimierung durch, um den optimalen Punkt zu markieren
options = optimset('Display', 'off');
optimal_values = fmincon(@(x) compute_dt(x), [0.02, 0.025], [], [], [], [], [0, 0], [0.04574, 0.04574], @(x)nonlcon(x), options);
optimal_dt = compute_dt(optimal_values);

% Berechne R und V mit den optimalen Werten
optimal_R = 2 * 0.45 / (1.225 * (0.12 * 0.4453^(2/3) + 0.86 * optimal_values(1)));
optimal_V = abs(nthroot(2 * 166.2 / (1.225 * (0.2403 * 0.4453^(2/3) + 0.6636 * optimal_values(1)) * cosd(16)), 3));


% Ausgabe der Ergebnisse
disp(['Optimaler Wert für S_p: ', num2str(optimal_values(1))]);
disp(['Optimaler Wert für S_r: ', num2str(optimal_values(2))]);
disp(['Minimiertes dt: ', num2str(optimal_dt)]);
disp(['Optimaler Wert für R: ', num2str(optimal_R)]);
disp(['Optimaler Wert für V: ', num2str(optimal_V)]);

% Funktion zur Berechnung von dt abhängig von S_p und S_r
function dt = compute_dt(x)
    S_p = x(1);
    S_r = x(2);

    R = 2 * 0.45 / (1.225 * (0.12 * 0.4453^(2/3) + 0.86 * S_p));
    V = abs(nthroot(2 * 166.2 / (1.225 * (0.25 * 0.4453^(2/3) + 0.6636 * S_p) * cosd(16)), 3));
    dt = R / V * 180 * pi / 180;
end
function [c, ceq] = nonlcon(x)
    S_p = x(1);
    S_r = x(2);
    R_min = 0.1; % Mindestwert für den Radius
    R_max = 13.0; % Maximalwert für den Radius

    R = 2 * 0.45 / (1.225 * (0.12 * 0.4453^(2/3) + 0.86 * S_p));

    % Nebenbedingungen
    c = [R - R_max; R_min - R; S_r + S_p - 0.04574]; % R muss zwischen R_min und R_max liegen, S_r + S_p muss gleich 0.04574 sein
    ceq = [];
end

