function sm = fn_set_MechParams(sm, perturb)
if(nargin < 2)
    perturb = 0;
end
% perturb - perturb the initial condition - perturb = 1 for first
% perturbation, perturb = 2 for second perturbation

%% Mechanical Parameters for ATRIAS Simulation
% All quantities are in SI Units

%% Environment
sm.g = 9.803;

%% Motors
% Saggital
sm.k_leg_spring = 3500;
sm.b_leg_spring = 0.075*sqrt(sm.k_leg_spring*0.5*.25^2);
sm.N = 50;
sm.mu_motor = 1.0;
sm.b_rotor = 0*60/(2*pi)*9.725e-5;
sm.i_rotor = 1.5e-3; 
sm.i_rotor_moments = [0 0 sm.i_rotor];
sm.i_motor_reflected = sm.i_rotor * sm.N^2;
sm.r_rotor = 0.001;
sm.m_rotor = 1e-5;
sm.r_virtual_gear = sm.N*sm.r_rotor;
sm.m_virtual_gear = 1e-5;
sm.i_virtual_gear = [0 0 1e-5];
sm.m_output_gear = 1e-5;
sm.i_output_gear = [0 0 1e-5];
sm.m_gearing = sm.m_rotor+sm.m_virtual_gear+sm.m_output_gear;
sm.Kt_leg = 0.121; % current constant
sm.Kb_leg = 0.121;
sm.R_leg = 0.052;
sm.L_leg = 0.102e-3;
sm.max_voltage_leg = 50;
sm.max_current_leg = 165*0.70;
sm.LEG_MTR_MAX_TORQUE = sm.max_current_leg * sm.Kt_leg * sm.N; % Nm
% Lateral
sm.r_hip_gearhead = 0.009525;
sm.r_hip_shaft = 0.542925;
sm.N_h = sm.r_hip_shaft / sm.r_hip_gearhead;
sm.mu_hmotor = 0.65;
sm.b_hmotor = 0;
sm.i_hrotor_moments = [0 0 3.193e-5];
sm.max_current_hip = 60;
sm.Kt_hip = 0.184;
sm.HIP_MTR_MAX_TORQUE = sm.max_current_hip*sm.Kt_hip*sm.N_h;



%% Motor Friction
% Vector order =  [RBack RFront LBack LFront]
% Specified after the gear reduction
sm.breakaway_friction = [25 25 25 25]*0.01; 
sm.coulomb_friction = [20 20 20 20]*0.01; 
sm.viscous_friction = [0.001 0.001 0.001 0.001];
sm.friction_linear_zone = [0.5 0.5 0.5 0.5].*pi/180; % velocity zone before full stiction is applied

%% Leg Geometry
sm.l_seg = 0.5;
sm.l_rad = 0.03;
% Thigh
sm.com_thigh = [0.0160 -0.1822 0];
sm.m_thigh = 0.69;
sm.i_thigh_moments = [0.02 0 0.02];
sm.i_thigh_products = [0 0 0];
% Lower Leg
sm.lowerleg_12 = 0.40;
sm.com_lowerleg = [0 -0.1513 0];
sm.m_lowerleg = 0.4438;
sm.i_lowerleg_moments = [0.01 0 0.01];
sm.i_lowerleg_products = [0 0 0];
% Four Bar Link
sm.com_fblink = [0 -0.1137 0];
sm.m_fblink = 0.46;
sm.i_fblink_moments = [0.01 0 0.01];
sm.i_fblink_products = [0 0 0];
% Shin
sm.shin_12 = 0.40;
sm.com_shin = [0.0215 -0.1503 0];
sm.m_shin = 0.75;
sm.i_shin_moments = [0.02 0 0.02];
sm.i_shin_products = [0 0 0];
% Mechanical Limits Model
sm.k_phi = 5000; % N/m
sm.v_phi_max = 1/60; % rad/s
sm.delta_phi_min = 25.5*pi/180;
sm.delta_phi_max = 160*pi/180;
sm.phi_shin_max = 3*pi/2 + 17.5*pi/180;
sm.phi_shin_min = pi - 45*pi/180;
sm.phi_thigh_max = pi + 45*pi/180;
sm.phi_thigh_min = pi/2 - 17.5*pi/180;
sm.phi_roll_max = pi + 16*pi/180;
sm.phi_roll_min = pi - 9*pi/180;
sm.phi_max_collision = -5*pi/180;

%% Hip Shaft
sm.leg_offset = 0.1831;
sm.hshaft_12 = [0 0 sm.leg_offset];
sm.com_hshaft = [0 0.029 sm.leg_offset-0.0149];
sm.m_hshaft = 18; % includes spring plates
sm.i_hshaft_moments = [0.29 0.10 0.27];
sm.i_hshaft_products = [0 0 0.1];

%% Boom setup
sm.pitch_mount_angle = 7.2824 * pi/180;
% boom tube
sm.l_boom = 1.8161;
sm.boom_12 = [0 0 -sm.l_boom];
sm.com_boom = [0 0 -1.33];
sm.m_boom = 4.49;
sm.i_boom_moments = [3.9 3.9 0];
sm.i_boom_products = [0 0 0];
% yaw shaft
sm.yshaft_12 = [0 0.9906 0];
sm.com_yshaft = [0 1 -0.15];
sm.m_yshaft = 12.89;
sm.i_yshaft_moments = [8 3.3 4.9];
sm.i_yshaft_products = [-3 0 0];
% boom base
sm.base_12 = [0 0.0286 0];
% counter weight
sm.counterweight_mass = 24.05;

%% Torso
sm.boom_mount_to_hip = 0.3176;
sm.boom_mount_to_center = 0.22339;
sm.boom_mount_to_center_diagonal = sm.boom_mount_to_center / cos(sm.pitch_mount_angle);
sm.torso_12 = [0 -sm.boom_mount_to_hip -sm.boom_mount_to_center];
sm.com_torso = [0 0.1824 -0.1577];
sm.m_torso = 21.9034;
sm.i_torso_moments = [1.5 1.5 2.2];
sm.i_torso_products = [0 0 0];
sm.winch_hook_offset = [0 0.375 -0.2225];
sm.d_hip_com_to_torso_com = sm.boom_mount_to_hip - sm.com_hshaft(2) + sm.com_torso(2);
sm.d_pelvis_to_IMU = [4.720 18.770 6.327]*2.54/100; 

%% Total masses and inertias
sm.num_legs = 2;
sm.m_total = sm.num_legs*(sm.m_hshaft + sm.m_thigh + sm.m_lowerleg + sm.m_shin + sm.m_fblink + 2*sm.m_gearing) + sm.m_torso;
sm.m_leg = sm.m_thigh + sm.m_lowerleg + sm.m_shin + sm.m_fblink;
sm.m_hip = sm.m_hshaft + sm.m_leg;
sm.i_motor = sm.i_motor_reflected + sm.m_leg/2*sm.l_seg^2;
sm.i_leg = sm.m_leg/2*sm.l_seg^2;
sm.d_com = (sm.d_hip_com_to_torso_com*sm.m_torso)/(sm.m_torso + sm.num_legs*sm.m_hshaft);
sm.boom_mount_to_com = sm.boom_mount_to_hip - sm.d_com;
sm.i_robot = sm.i_torso_moments(3) + sm.m_torso*(sm.d_hip_com_to_torso_com - sm.d_com)^2 + sm.num_legs*sm.m_hshaft*sm.d_com^2 + sm.num_legs*sm.i_hshaft_moments(1);

%% Winch
sm.winch_gravity_comp = sm.m_total*sm.g*sm.l_boom;
sm.kp_winch = sm.winch_gravity_comp/(3*pi/180);
sm.td_winch = 50e-3;
sm.kd_winch = sm.td_winch*sm.kp_winch;

%% Initial parameters
% -------------------
sm.initial_l_leg_length = 0.86; %0.96;
sm.initial_l_leg_angle = pi/2 + pi/2;
sm.initial_r_leg_length = 0.7;
sm.initial_r_leg_angle = pi/2 + pi/2;
sm.thigh_initial_r = sm.initial_r_leg_angle - acos(sm.initial_r_leg_length);
sm.shin_initial_r = sm.initial_r_leg_angle + acos(sm.initial_r_leg_length);
sm.thigh_initial_l = sm.initial_l_leg_angle - acos(sm.initial_l_leg_length);
sm.shin_initial_l = sm.initial_l_leg_angle + acos(sm.initial_l_leg_length);
sm.v_thigh_initial_r =  0;
sm.v_shin_initial_r =  0;
sm.v_thigh_initial_l =  0;
sm.v_shin_initial_l =  0;
sm.initial_boom_yaw = 0;
sm.initial_boom_pitch = 0;
sm.initial_com_height = 0.94; %0.94; %1.2;%1.0
sm.initial_torso_height = sm.initial_com_height + sm.boom_mount_to_com;
sm.initial_boom_roll = asin((sm.initial_com_height + sm.boom_mount_to_com + sm.boom_mount_to_center*tan(sm.pitch_mount_angle) - (sm.base_12(2) + sm.yshaft_12(2))) / (sm.l_boom + sm.boom_mount_to_center_diagonal));
sm.initial_winch_angle = -12 * pi/180;%-12 * pi/180;
sm.initial_hip_roll = 0*pi/180;
sm.initial_dx = 0;
sm.initial_boom_yaw_velocity = sm.initial_dx / ((sm.l_boom+sm.boom_mount_to_center_diagonal)*cos(sm.initial_boom_roll));

%initial deflection to load the left leg
sm.initial_deflection_l = 0;%0.5*sm.m_total*sm.g*sin((sm.shin_initial_l-sm.thigh_initial_l)*0.5)/sm.k_leg_spring;

%% Ground Contact Model
% -------------------
% softening constant
sm.soften_const = 1.0;

% vertical ground interaction stiffness
sm.k_gy = sm.m_total*sm.g/0.005; % [N/m]

% maximum vertical ground relaxation speed
sm.v_gy_max = 0.03; %[m/s]

% horizontal ground interaction stiffness
sm.k_gx = sm.m_total*sm.g/0.01; % [N/m]

% maximum horizontal ground relaxation speed
sm.v_gx_max = 0.03; %[m/s]

% stiction coefficient
sm.mu_stick = 0.9;

% sliding coefficient
sm.mu_slide = 0.8;

% slip-stic transition velocity
sm.vLimit = 0.01; %[m/s]

% ground height
sm.floor_tile_width = 30*pi/180;
sm.max_floor_tile_height =0*0.96*2.0; %10cm
if perturb
sm.max_floor_tile_height =10/100; %10cm
end
sm.tile_cross_section = [1.5 0.03; 1.5 0; 3 0; 3 0.03];
sm.tile_delta = [0 1; 0 0; 0 0; 0 1];
% rng(156); % random seed
tile_scale = rand(2*pi/sm.floor_tile_width,1); % random heights
tile_scale(1:2) = 0;
tile_scale(end) = 0;
sm.tile_scale = tile_scale;
sm.tile_scale = sm.max_floor_tile_height * sm.tile_scale / max(abs(sm.tile_scale));
%tile_scale = [0 8.2 6.9 12.6 0 20 20 0 0 12.1 15.5 6.3]'./100; % manual heights
sm.max_floor_tile_height = max(sm.tile_scale);
% tile_scale(1) = 0;
sm.cmap = jet(128);
sm.tile_color = sm.cmap(floor((sm.tile_scale/(sm.max_floor_tile_height+eps))*(length(sm.cmap)-1))+1,:);