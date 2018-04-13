function [mc, idc, rc] = set_control_params_raibert_id(sm, sample_time, disturb)
if(nargin < 3)
    disturb = 0;
end
% mc - model control parameters - sensors, sea, filters
% idc - inverse dynamics parameters - lqr, etc
% rc - raibert control parameters
%% Set control parameters for sensors, filters, SEA
mc.g = sm.g;
mc.k_leg_spring = sm.k_leg_spring;

% Sensors
mc.fcut_filter1 = 400;
mc.fcut_acceleration = 600;

% SEA Control
mc.kp_sea_torque = 10;
mc.kd_sea_torque = 0.9*2*sqrt(mc.kp_sea_torque/3500*3.75);
mc.k_ff = 0.75;

% optimal motor torque distribution
mc.leg_motor_saturation = sm.LEG_MTR_MAX_TORQUE;
mc.lambda1 = 1; % hip torque 
mc.lambda2 = 1;  % leg force

% Hip roll control
mc.m_leg_motor = 18*(1 + disturb*(rand-0.5)); % includes spring plates
mc.m_hshaft = mc.m_leg_motor;
mc.lateral_offset = 0.1831;
mc.leg_offset = mc.lateral_offset;
mc.hip_gravity_compensation = mc.m_hshaft*mc.g*mc.leg_offset; % hip_mass * gravity * lever_distance
mc.max_torque_lateral = 4 * mc.hip_gravity_compensation; 

mc.HIP_MTR_MAX_CONT_CURRENT = 8.2; %Maximum continuous hip motor current
mc.HIP_MOTOR_CONSTANT = 0.184;
mc.HIP_MTR_GEAR_RATIO = 57;
mc.HIP_MTR_MAX_CONT_TORQUE = mc.HIP_MTR_MAX_CONT_CURRENT*mc.HIP_MOTOR_CONSTANT*mc.HIP_MTR_GEAR_RATIO;

mc.max_primary_lateral_stance_torque = 2*mc.HIP_MTR_MAX_CONT_TORQUE;
mc.kp_lateral = 8000;
mc.kd_lateral = 80;

% Leg Length Control
mc.kp_leg_length = 1000;
mc.td_leg_length = 40e-3;
mc.kd_leg_length = mc.td_leg_length*mc.kp_leg_length;

% Leg Angle Control
mc.kp_leg_angle = 900;
mc.kd_leg_angle = 90;

% Torso Balance Control
% desired_torso_angle = 0*pi/180;
mc.kp_torso = 1000;
mc.td_torso = 200e-3;
mc.kd_torso = mc.td_torso*mc.kp_torso;

% contact detection
mc.fcut_contact = 600;
mc.contact_threshold = 0.3;
mc.loaded_threshold = 0.50;


% Velocity-based swing trajectory following
mc.kp_velocity_swing_trajectory = 10;

%% Raibert Hopper Parameters

% Spring parameters
rc.l0_virtual = 0.90;
rc.kp_virtual = 13000;
rc.kd_virtual = 1e-3*rc.kp_virtual;
rc.kp_flight_max = 20000;
rc.kp_flight_min = 5000;
rc.kd_flight = 250;
rc.l_min = 0.4;
rc.l_retract = 0.90;
rc.z_apex = 0.15;

% placement control
rc.alpha_min = 60*pi/180;
rc.x_limit = rc.l0_virtual*cos(rc.alpha_min);
rc.k_placement = 0.1;
rc.t_alpha = 1.0;
rc.T_step = 0.3;
rc.desired_xd = 1.0;
rc.desired_theta = 0.0;
rc.desired_z = 0.95;
rc.kpz = 3000;
rc.kdz = 300;

rc.kpx = 0;
rc.kdx = 0;
rc.kpt = -1000;
rc.kdt = -150;

rc.gain_ff = 1.0;
rc.Cd = 0;

rc.dx_avg_samples = 0.03/sample_time;

% State Transition triggers
rc.l_max = 0.95; %m (threshold to stop thrusting)
rc.t_flight = 0.001; %[s] (escape twilight period)

% thrust control
rc.F_thrust =   300; %[N] (leg thrust force)
rc.dt_thrust =  0.080; %[s] (thrust period)

rc.manual_leg_selection = 1;

desired_standing_z = 0.9;
y_star = 0.15;
rc.X_Star_left = [0; 0; 0];
rc.Y_Star_left = [y_star; 0; 0];
rc.Z_Star_left =  [desired_standing_z; 0; 0];

rc.X_Star_right = [0; 0; 0];
rc.Y_Star_right = [-y_star; 0; 0];
rc.Z_Star_right =  [desired_standing_z/2; 0; 0];
%%
mc.lpf_damping = sqrt(2)/2; % Butterworth filter damping
mc.fcut_filter1 = 20*(2*pi); % Cutoff frequency for velocities
mc.B1_lpf_filter1 = -2*exp(-mc.lpf_damping*mc.fcut_filter1*sample_time)*cos(mc.fcut_filter1*sample_time*sqrt(1-mc.lpf_damping^2));
mc.B2_lpf_filter1 = exp(-2*mc.lpf_damping*mc.fcut_filter1*sample_time);
mc.A_lpf_filter1 = 1 + mc.B1_lpf_filter1 + mc.B2_lpf_filter1;


%% Leg Limits for commands
mc.min_front_bar_angle = 98.0*pi/180; % ~83 degrees
mc.max_front_bar_angle = 194.8*pi/180; % ~218 degrees
mc.min_back_bar_angle = 165.0*pi/180; % ~144 degrees
mc.max_back_bar_angle = 262.4*pi/180; % ~279 degrees
mc.max_leg_length = 0.95;
mc.min_leg_length = 0.5;

mc.kp_swing_trajectory = 2500;
mc.kd_swing_trajectory = 145;

%% Boom
mc.m_boom = (4.49 + 3)*(1 + disturb*(rand-0.5)); % kg
mc.l_boom = 1.767; % m, length of boom
mc.h_boom = 0.99; % m, height of boom center

%% Lateral Motion
mc.lateral_motor_efficiency = 0.65;
mc.lateral_offset = 0.1831; % m
mc.r_hip_gearhead = 0.009525;
mc.r_hip_shaft = 0.542925;
mc.N_hip = mc.r_hip_shaft / mc.r_hip_gearhead;

mc.k_sea_low = [3343; 3825; 3476; 3905]; 
mc.k_sea_high = [4255.0; 4255.0; 4372.1; 4322.4];
mc.k_sea = mc.k_sea_low;
mc.k_leg_spring = mean(mc.k_sea_low);

mc.LEG_MTR_MAX_TORQUE = sm.LEG_MTR_MAX_TORQUE; 
mc.LEG_MTR_GEAR_RATIO = sm.N;
% mc.m_leg_motor = 18; % kg
mc.com_leg_motor = [0, mc.lateral_offset-0.0149, 0.029]; % m, Coordinates from pelvis for left leg. Flip y coordinate for right leg.
mc.i_leg_motor = [0.29, 0.27,  0.10];
mc.j_rotor = 1.5e-3;
mc.j_motor = mc.j_rotor*mc.LEG_MTR_GEAR_RATIO^2;
mc.j_leg = 2*0.147;
mc.j_segments = [0.145 0.155 0.145 0.155];
mc.leg_motor = 18; % kg
mc.com_leg_motor = [0, mc.lateral_offset-0.0149, 0.029]; % m, Coordinates from pelvis for left leg. Flip y coordinate for right leg.
mc.i_leg_motor = [0.29, 0.27,  0.10];
mc.j_rotor = 1.5e-3;
mc.j_motor = mc.j_rotor*mc.LEG_MTR_GEAR_RATIO^2;
mc.j_leg = 2*0.147;
mc.j_segments = [0.145 0.155 0.145 0.155];

%% Legs
mc.m_leg = 2.3438*(1 + disturb*(rand-0.5)); % kg
mc.l_seg = 0.5; % m
mc.l_seg_short = 0.4; % m
mc.d_ankle_to_foot = 0.0245; % m
mc.l_seg_lower = 0.09 + mc.d_ankle_to_foot; % m

%% Full Robot
mc.m_total_real = 62.591*(1 + disturb*(rand-0.5)); % kg
mc.m_total_real = mc.m_total_real + mc.m_boom/2;
mc.d_vertical_mount_to_hip = 0.3176; % m

%% Torso
mc.m_torso = mc.m_total_real - 2*(mc.m_leg_motor+mc.m_leg);
mc.m_body = mc.m_total_real - 2*mc.m_leg;
mc.com_torso = [0, 0.0, 0.50];%[0.01, 0, 0.2885]; % m, Coordinates from pelvis.
%if robot_is_attached_to_boom, com_torso(2) = -0.1; end
mc.i_torso = [1.5, 2.2, 1.5]*(1 + disturb*(rand-0.5));
mc.pelvis_to_imu = [0.11988, -0.16071+0.07366/2, 0.47675];
%% Extra boom quantities
mc.l_yaw_center_to_pelvis = 1.978; % m, distance from center of boom circle to unpitched pelvis
mc.boom_pelvis_offset_angle = 0.174439; % rad, angle offset between above line and boom roll
mc.boom_mount_angle = 7.2824 * pi/180; % rad, angle offset due to boom mount
mc.d_horizontal_mount_offset = 0.22339; % m, horizontal distance from center of boom attachment to frontal centerline
mc.d_vertical_mount_offset = mc.d_horizontal_mount_offset * tan(mc.boom_mount_angle); % m, vertical distance from center of boom attachment to frontal centerline
mc.l_boom_projected = mc.l_boom + norm([mc.d_horizontal_mount_offset, mc.d_vertical_mount_offset]) ; % distance of boom projected to frontal centerline
mc.d_vertical_boom_hip_offset = mc.d_vertical_mount_to_hip + mc.d_vertical_mount_offset; % distance from (boom attachment projected onto frontal center) to (hip rotation point)
mc.boom_mount_to_center_diagonal =  norm([mc.d_horizontal_mount_offset, mc.d_vertical_mount_offset]) ;  

%% Approximate (static) total system COM distances and inertia

mc.d_vertical_com = norm(mc.com_torso([1 3]))*mc.m_torso/mc.m_total_real; % m, Distance above saggital rotation point. 
mc.d_vertical_mount_to_com = mc.d_vertical_mount_to_hip - mc.d_vertical_com; % m
mc.d_vertical_boom_com_offset = mc.d_vertical_mount_offset + mc.d_vertical_mount_to_com; % distance from (boom attachment projected onto frontal center) to (system center of mass)
mc.d_vertical_com_to_trunk_com = norm(mc.com_torso([1 3]))-mc.d_vertical_com;
mc.i_robot = 2*(mc.m_leg_motor+mc.m_leg)*mc.d_vertical_com^2 + mc.m_torso*mc.d_vertical_com_to_trunk_com^2 + mc.i_torso(2) + 2*mc.i_leg_motor(2); % kg m^2, y principal component

%% SEA filters
mc.fcut_dtau = 80*(2*pi); % Hz, Low pass filter cutoff frequency when calculating acceleration from velocity
mc.lpf_damping = sqrt(2)/2; % butterworth damping ratio
mc.B1_lpf_dtau = -2*exp(-mc.lpf_damping*mc.fcut_dtau*sample_time)*cos(mc.fcut_dtau*sample_time*sqrt(1-mc.lpf_damping^2));
mc.B2_lpf_dtau = exp(-2*mc.lpf_damping*mc.fcut_dtau*sample_time);
mc.A_lpf_dtau = 1 + mc.B1_lpf_dtau + mc.B2_lpf_dtau;

%% -- Low pass filter for IMU measurements
mc.fcut_imu = 250*(2*pi);
mc.lpf_damping = sqrt(2)/2;
mc.B1_lpf_imu = -2*exp(-mc.lpf_damping*mc.fcut_imu*sample_time)*cos(mc.fcut_imu*sample_time*sqrt(1-mc.lpf_damping^2));
mc.B2_lpf_imu = exp(-2*mc.lpf_damping*mc.fcut_imu*sample_time);
mc.A_lpf_imu = 1 + mc.B1_lpf_imu + mc.B2_lpf_imu;
mc.fcut_accelerometer = 80*(2*pi);
mc.lpf_damping = sqrt(2)/2;
mc.B1_lpf_accelerometer = -2*exp(-mc.lpf_damping*mc.fcut_accelerometer*sample_time)*cos(mc.fcut_accelerometer*sample_time*sqrt(1-mc.lpf_damping^2));
mc.B2_lpf_accelerometer = exp(-2*mc.lpf_damping*mc.fcut_accelerometer*sample_time);
mc.A_lpf_accelerometer = 1 + mc.B1_lpf_accelerometer + mc.B2_lpf_accelerometer;

%% Lateral Kalman filter

mc.fcut_lateral_angle = 80*(2*pi);
mc.lpf_damping = sqrt(2)/2;
mc.B1_lpf_lateral_angle = -2*exp(-mc.lpf_damping*mc.fcut_lateral_angle*sample_time)*cos(mc.fcut_lateral_angle*sample_time*sqrt(1-mc.lpf_damping^2));
mc.B2_lpf_lateral_angle = exp(-2*mc.lpf_damping*mc.fcut_lateral_angle*sample_time);
mc.A_lpf_lateral_angle = 1 + mc.B1_lpf_lateral_angle + mc.B2_lpf_lateral_angle;

%% Kalman Filtering parameters - Center of Mass
% Coordinate system is z-up and in the yawed frame.
% x filter states: [dw_x]
% s_k+1 = s_k + dt/m*F
% when in stance, measure dx_w ~ dx_f with confidence based on loading
mc.A_kalman_transverse = [1 sample_time 0; 0 1 sample_time; 0 0 1];
mc.B_kalman_transverse = [0; 0; 1/mc.m_total_real];
mc.C_kalman_transverse = [1 0 0; 1 0 0; 0 0 1];
mc.G_kalman_transverse = mc.B_kalman_transverse;
mc.A_kalman_no_position = [1 sample_time; 0 1];
mc.B_kalman_no_position = [0; 1/mc.m_total_real];
mc.C_kalman_no_position = [1 0; 1 0; 0 1];
mc.G_kalman_no_position = mc.B_kalman_no_position;
mc.A_kalman_vertical = [1 sample_time 0; 0 1 sample_time; 0 0 1];
mc.B_kalman_vertical = [0; 0; 1/mc.m_total_real];
mc.C_kalman_vertical = [1 0 0; 1 0 0; 0 0 1];
mc.G_kalman_vertical = mc.B_kalman_vertical;
% Covariances
mc.Q_kalman_GRF = 20^2; % N, error in GRF estimates
mc.Q_kalman_GRFz_difference = 2*(50)^2 / sample_time;
mc.Q_kalman_GRFx_difference = 2*(50)^2 / sample_time;
mc.Q_kalman_GRFy_difference = 2*(50)^2 / sample_time;
mc.R_kalman_accelerometer = 1.00^2;  % m/s^2
mc.R_kalman_foot_stance = 0.01^2 + 0.005^2;  % m, geometry error + ground error
mc.R_kalman_foot_slip = 0.01^2 + 1^2;
mc.R_kalman_dfoot_stance = 0.25^2; % m/s
mc.R_kalman_dfoot_slip = 5^2; % m/s
mc.kalman_stance_slip_threshold = 0.1;
% Initial estimates
mc.P0_kalman_transverse = diag([0.015^2, 0.05^2, 0.125^2]); % [m, m/s, m/s^2], initial state error covariance
mc.P0_kalman_no_position = diag([0.05^2, 0.125^2]); % [m/s, m/s^2], initial state error covariance
mc.P0_kalman_vertical = diag([0.015^2, 0.05^2, 0.125^2]); % [m, m/s, m/s^2], initial state error covariance

%% Inverse Dynamics stuff
idc.desired_torso_pitch = 0*pi/180;
idc.desired_torso_angle = idc.desired_torso_pitch;
idc.desired_z = 0.95;



%% PD control on swinging hip
rc.kp_lateral_beta2 = 1000;
mc.lateral_hip_inertia_about_pelvis = (mc.m_leg+mc.m_leg_motor)*mc.lateral_offset^2;
rc.kd_lateral_beta2 = 0.8*2*sqrt(rc.kp_lateral_beta2*mc.lateral_hip_inertia_about_pelvis);
rc.k_stance_secondary_roll = [rc.kp_lateral_beta2; rc.kd_lateral_beta2];
rc.k_lateral_swing_on_boom = rc.k_stance_secondary_roll;

end