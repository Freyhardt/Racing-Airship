%% Author: Jan S. Freyhardt \n,
%% Bachelor Thesis \n,
%% Evelope Calculation Positions Calculating Tool
%% email: jan.freyhardt@gmail.com\n

clc
clear
close all


g = 9.81; 

% Atmospheric Conditions
T_0 = 273.15; % Temperature in Kelvin [K]
T_ISA = 288.15; % 15° ISA standard Temperature in Kelvin 
R_Air = 287.1; % Specific gas constant for dry air [J/(kg*K)]
p_0 = 101325; % Pressure in [Pa]
Rho_0 = p_0/(R_Air*T_0);

L = 3.008310292; % Envelope length in millimeters [mm]
D_max = 0.546965508; % Envelope diameter in millimeters [mm]
EnvelopeVolume = 0.44532; % in cubic millimeters [m^3]
M_Envelope = 0.12047; % Mass Envelope in [kg]

% Envelope modeling parameters
a_1 = 1.4;
a_2 = -1.75034722222225;
a_3 = -3.26238425925885;
a_4 = 11.2749131944435;
a_5 = -12.7612847222211;
a_6 = 5.0991030092589;

% Normalize x values between 0 and 1 for positive values
x_normalized_positive = linspace(0, 1, 1000);

% Envelope Equation
EnvelopeEquation_positive = D_max * sqrt(a_1*x_normalized_positive + a_2*x_normalized_positive.^2 + a_3*x_normalized_positive.^3 + a_4*x_normalized_positive.^4 + a_5*x_normalized_positive.^5 + a_6*x_normalized_positive.^6);

% Envelope Figure
figure;
plot(x_normalized_positive * L, real(EnvelopeEquation_positive), 'b', x_normalized_positive * L, -real(EnvelopeEquation_positive), 'r')
xlabel('Z Axis [mm]') % Provide appropriate x-axis label
ylabel('X Axis [mm] ') % Provide appropriate y-axis label
title('Envelope Model') % Provide a title for the plot

% Set axis equal
axis equal;

% Previous lines remain unchanged

% Interpolation of the envelope data for positive x values
x_normalized_positive = linspace(0, 1, 1000);
EnvelopeEquation_positive = D_max * sqrt(a_1*x_normalized_positive + a_2*x_normalized_positive.^2 + a_3*x_normalized_positive.^3 + a_4*x_normalized_positive.^4 + a_5*x_normalized_positive.^5 + a_6*x_normalized_positive.^6);

% Create interpolation function for positive x values
envelopeInterpolation_positive = @(x) interp1(x_normalized_positive * L, real(EnvelopeEquation_positive), x, 'linear', 'extrap');

% User input for x value
userInputX_positive = input('Enter the positive x value: ');

% Obtain y value through interpolation
userOutputY_positive = envelopeInterpolation_positive(userInputX_positive);

% Display the result
disp(['For x = ', num2str(userInputX_positive), ', the corresponding positive y value is: ', num2str(userOutputY_positive)]);

% Interpolation of the envelope data for negative x values
x_normalized_negative = -linspace(0, 1, 1000);
EnvelopeEquation_negative = -D_max * sqrt(a_1*x_normalized_negative + a_2*x_normalized_negative.^2 + a_3*x_normalized_negative.^3 + a_4*x_normalized_negative.^4 + a_5*x_normalized_negative.^5 + a_6*x_normalized_negative.^6);

% Create interpolation function for negative x values
envelopeInterpolation_negative = @(x) interp1(x_normalized_negative * L, real(EnvelopeEquation_negative), x, 'linear', 'extrap');

% User input for x value
userInputX_negative = input('Enter the negative x value: ');

% Obtain y value through interpolation
userOutputY_negative = envelopeInterpolation_negative(userInputX_negative);

% Display the result
disp(['For x = ', num2str(userInputX_negative), ', the corresponding negative y value is: ', num2str(userOutputY_negative)]);


% Cross-sectional Area
ProfileArea = integral(@(x) (sqrt(a_1*x + a_2*x.^2 + a_3*x.^3 + a_4*x.^4 + a_5*x.^5 + a_6*x.^6)), 0, L);

% Volume Calculation
% Volume = integral(@(x) D_max * sqrt(a_1*x_normalized_positive + a_2*x_normalized_positive.^2 + a_3*x_normalized_positive.^3 + a_4*x_normalized_positive.^4 + a_5*x_normalized_positive.^5 + a_6*x_normalized_positive.^6),0,L)^2* pi;

% Volume centroid
X_centroid = integral(@(x) (sqrt(a_1*x + a_2*x.^2 + a_3*x.^3 + a_4*x.^4 + a_5*x.^5 + a_6*x.^6) .* x), 0, L) / (2*ProfileArea);

%disp(['Volume: ', num2str(Volume)]);
disp(['Volume centroid: ', num2str(X_centroid)]);


% Parameters Lifting gas
p_diff = 100; % Temperature differential

R_Gas = 2077.1;
RhoGas_0 = 0.17; 


% Lift
H = 400; % Height in [m]
p_H = p_0*exp(-(g*H)/(R_Air*T_ISA)); % Pressure at height in [Pa]
p_Gas = p_H+p_diff; 


M_Lift = (p_H/(R_Air*T_ISA)-p_Gas/(R_Gas*T_ISA)) *EnvelopeVolume; % Lifting mass in [kg]
F_Lift = M_Lift*g; % Envelope lifting force in [N]
M_LiftTotal = M_Lift-M_Envelope; % Total Lifting force of the Envelope

