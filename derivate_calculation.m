%% Author: Jan S. Freyhardt \n,
%% Bachelor Thesis \n,
%% Derivate Calculation
%% email: jan.freyhardt@gmail.com\n

clc
clear
close all



g = 9.81; 

%Atmosparic Conditions
T_0 = 273.15; % Temperature in Kalvin [K]
T_ISA = 288.15; % 15° ISA satandart Temperature in Kalvin 
R_Air = 287.1; % Specific gas constant for dry air [J/(kg*K)]
p_0 = 101325; % Pressure in [Pa]
Rho_0 = p_0/(R_Air*T_0);

L = 3.008310292; % Envelope length in milimeters [mm]
D_max = 0.546965508; % Envelope dimeter in milimeters [mm]
EnvelopeVolume = 0.44532; % in cubicmilimeters [m^3]
M_Envelope = 0.12047; % Mass Envelope in [kg]

% Envelope modelation parameters
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
xlabel('Z Axis [m]') % Provide appropriate x-axis label
ylabel('X Axis [m] ') % Provide appropriate y-axis label
title('Envelope Model') % Provide a title for the plot

% Set axis equal
axis equal;

% Querschnittsfläche (Profile Cross-sectional Area)
ProfileArea = integral(@(x) (sqrt(a_1*x + a_2*x.^2 + a_3*x.^3 + a_4*x.^4 + a_5*x.^5 + a_6*x.^6)), 0, L);

% Volumenberechnung
% Volumen = integral(@(x) D_max * sqrt(a_1*x_normalized_positive + a_2*x_normalized_positive.^2 + a_3*x_normalized_positive.^3 + a_4*x_normalized_positive.^4 + a_5*x_normalized_positive.^5 + a_6*x_normalized_positive.^6),0,L)^2* pi;

% Volumenschwerpunkt
X_schwerpunkt = integral(@(x) (sqrt(a_1*x + a_2*x.^2 + a_3*x.^3 + a_4*x.^4 + a_5*x.^5 + a_6*x.^6) .* x), 0, L) / (2*ProfileArea);

%disp(['Volumen: ', num2str(Volumen)]);
disp(['Volumenschwerpunkt: ', num2str(X_schwerpunkt)]);


%Parameters Lifftinggas
p_diff = 100; % Temperature differential

R_Gas = 2077.1;
RhoGas_0 = 0.17; 


%Lift
H = 400; % Hight in [m]
p_H = p_0*exp(-(g*H)/(R_Air*T_ISA)); % Pressure at hight in [Pa]
p_Gas = p_H+p_diff; 


M_Lift = (p_H/(R_Air*T_ISA)-p_Gas/(R_Gas*T_ISA)) *EnvelopeVolume; % Lifting mass in [kg]
F_Lift = M_Lift*g; % Envelope lifting force in [N]
M_LiftTotal = M_Lift-M_Envelope; % Total Liftigforce of the Envelope



%Finns
Munkfaktor = 0;
% Electronics

% Derivatives 

m = 20;
% Definiere die Variablen und Funktionen
x = sym('x');
theta = 0; % Der Nickwinkel ist immer gleich null im stationären Flug
D_max = 0.546965508;
Betha = 6 * (pi/180);
a = L/2; % Halbachse
b = sqrt(3/(4*pi*EnvelopeVolume));  % Radius
k_1 = 0.08; % k_1 Munkfaktor
k_2 = 0.86; % k_2 Munkfaktor
k_3 = 0.86; % k_3 Munkfaktor
K_1 = k_1 + 1; 
K_2 = k_2 + 1;
K_3 = (a^2+b^2)/(a^2-b^2);
Eta_k = 1; % fluid koefficient
Eta_f = 1;
Velocity = 21;
Lambda_f = 3;
c_l_alpha0 = (2 * pi * Lambda_f)/(sqrt(Lambda_f^2 + 4 +2));
S_f = 1;
x_f = 1.5;
l_R = 2;
r = D_max * sqrt(a_1*x + a_2*x^2 + a_3*x^3 + a_4*x^4 + a_5*x^5 + a_6*x^6);
A = pi*r^2;


% Definiere die zu integrierenden Funktionen
f = 1/4 * K_1 * (cos(theta)^2 * diff(A, x) - 2 * diff(A, x));
I_uw_z = 1/2 * K_1 * K_2 * cos(theta)^2  *diff(A, x);
I_uq_z = A + K_1 * (K_3-1) * cos(theta)^2 * A + 1/2 * K_1 * (K_3+1) * x * cos(theta)^2 * diff(A, x);
I_uw_m = 1/2 * K_1 * K_2 *(r * sin(theta)*cos(theta) * diff(A, x) + x * cos(theta)^2 * diff(A, x));
I_uq_m = 1/2 * r^2 * diff(A, x) + x * A + 1/2 * K_1 *(K_3-1) * (r^2 * cos(theta)^2 * diff(A, x)+ 2 * x * cos(theta)^2*A)+1/2*K_1*(K_3+1)*(x*r*sin(theta)*cos(theta)*diff(A, x)+x^2*cos(theta)^2*diff(A, x));

% Wandle die symbolischen Funktionen in numerische Funktionen um
numeric_function = matlabFunction(f);
numeric_function1 = matlabFunction(I_uw_z);
numeric_function2 = matlabFunction(I_uq_z);
numeric_function3 = matlabFunction(I_uw_m);
numeric_function4 = matlabFunction(I_uq_m);

% Berechne die Integrale numerisch
Ixpot_numeric = integral(numeric_function, 0, 3);
Ixpot_numeric1 = integral(numeric_function1, 0, 3);
Ixpot_numeric2 = integral(numeric_function2, 0, 3);
Ixpot_numeric3 = integral(numeric_function3, 0, 3);
Ixpot_numeric4 = integral(numeric_function4, 0, 3);

% Zeige die numerischen Ergebnisse an
disp('Das numerische Ergebnis des Integrals ist:');
disp(Ixpot_numeric);
disp('Das numerische Ergebnis von I_uw_z ist:');
disp(Ixpot_numeric1);
disp('Das numerische Ergebnis von I_uq_z ist:');
disp(Ixpot_numeric2);
disp('Das numerische Ergebnis von I_uw_m ist:');
disp(Ixpot_numeric3);
disp('Das numerische Ergebnis von I_uq_m ist:');
disp(Ixpot_numeric4);

c_y_betha = (2*Eta_k*Ixpot_numeric1-Eta_f*c_l_alpha0*S_f)/EnvelopeVolume^(2/3);
disp(c_y_betha);
c_y_r     = (2*Eta_k*Ixpot_numeric2-Eta_f*c_l_alpha0*S_f*x_f)/(EnvelopeVolume^(2/3)*l_R);
disp(c_y_r);
c_n_betha = (2*Eta_k*Ixpot_numeric3-Eta_f*c_l_alpha0*S_f*x_f)/(EnvelopeVolume^(2/3)*l_R);
disp(c_n_betha);
c_n_r     = (2*Eta_k*Ixpot_numeric3-Eta_f*c_l_alpha0*S_f*x_f^2)/(EnvelopeVolume^(2/3)*l_R^2);
disp(c_n_r);

R = ((2*m*x_f)/(Rho_0*EnvelopeVolume*l_R)-c_y_r*x_f+c_n_r*l_R)/(c_y_betha*(x_f/l_R)-c_n_betha)*(1/Betha);









