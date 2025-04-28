% Load parameters

M = 2200; % Total car mass (assuming motor mass is included)
g = 9.81; % Gravitational acceleration
mu_r = 0.01; % Rolling resistance coefficient
c_w = 0.2; % Wind resistance coefficient
rho = 1.2; % Air density
r_w = 0.25; % Wheel radius in meters
i_g = 10; % Gear ratio (Motor : Wheel)
A = 2.5; % Frontal area

% Inertia
J_g = 0.1; % Gear inertia ref motor shaft
J_m = 0.06; % Inertia of motor shaft
J_car = (((M/4)*r_w^2)/i_g^2); % Inertia of car
J = J_car + J_m + J_g; % Total inertia referred to motor shaft

% TTL250CV Motor parameters

% Armature
I_an = 160;
R_a = 0.038;
U_an = 190;
T_a = 0.200;
L_a = T_a*R_a;

% Field
I_fn = 9.5;
R_f = 7.9;
U_fn = I_fn * R_f;
T_f = 0.055;
L_f = T_f*R_f;

N_n = 4300; % Nominal rpm

% Scaled motor parameters
r_a = R_a*I_an/(U_an-R_a*I_an);
l_a = L_a*I_an/(U_an-R_a*I_an);

% Other time constants

T_samp = 0.000125;
T_f_f = 0;
T_f_a = 0;

% Controller parametes

U_dc = 350; % DC link voltage
K_pf = 320; % Field controller gain
K_pa = 3.1; % Armature controller gain
K_ps = 20; % Speed controller gain
