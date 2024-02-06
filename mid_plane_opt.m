%% Author: Jan S. Freyhardt \n,
%% Bachelor Thesis \n,
%% Mid Plane Optimisation
%% email: jan.freyhardt@gmail.com\n

clc
clear
close all

L = 3.008310292; % Envelope length in millimeters [mm]
D_max = 0.546965508; % Envelope diameter in millimeters [mm]
EnvelopeVolume = 0.44532; % in cubic millimeters [m^3]
M_Envelope = 0.12047; % Mass Envelope in [kg]

a_1 = 1.4;
a_2 = -1.75034722222225;
a_3 = -3.26238425925885;
a_4 = 11.2749131944435;
a_5 = -12.7612847222211;
a_6 = 5.0991030092589;

x_normalized_positive =  0.518; % You need to define the value of x_normalized_positive

% Initial guess for z_MP
z_MP = 0.2735;

% Iterative optimization
for iter = 1:10
    z = D_max * sqrt(a_1 * x_normalized_positive + a_2 * x_normalized_positive.^2 + a_3 * x_normalized_positive.^3 + a_4 * x_normalized_positive.^4 + a_5 * x_normalized_positive.^5 + a_6 * x_normalized_positive.^6);
    z_CG = (120.47 * 0 + 2 * 83.4 * 0 + 3 * 15 * 0 + 89.7 * z + 12 * z_MP + 8 * 0) / (120.47 + 2 * 83.4 + 3 * 15 + 89.7 + 12 + 8);
    S_MP = (2 * 0.45 * z_CG) / (1.225 * 0.86 * 9 * z_MP);
    
    % Update z_MP using the new S_MP
    z_MP = 0.2735 + (sqrt(2 * S_MP))/2;
    
    % Check for convergence (you may need to adjust the tolerance)
    if abs(S_MP - z_MP) < 1e-6
        break;
    end
end

disp(['Optimal S_MP: ' num2str(S_MP)]);
disp(['Optimal x: ' num2str(x_normalized_positive)]);
disp(['Corresponding z_MP: ' num2str(z_MP)]);
disp(['Corresponding z_CG: ' num2str(z_CG)]);

