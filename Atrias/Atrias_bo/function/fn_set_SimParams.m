function [sm, sample_freq, sample_time] = fn_set_SimParams(sm)

sample_freq = 1000;
sample_time = 1/sample_freq;
update_freq = sample_freq;
%% Timings
sm.t_end = 30;  
sm.output_sample_time = 1/120;

% Sensors
sm.lpf_damping = sqrt(2)/2; % Butterworth filter damping
sm.fcut_velocity = 100*(2*pi); % Cutoff frequency for velocities
sm.B1_lpf_velocity = -2*exp(-sm.lpf_damping*sm.fcut_velocity*sample_time)*cos(sm.fcut_velocity*sample_time*sqrt(1-sm.lpf_damping^2));
sm.B2_lpf_velocity = exp(-2*sm.lpf_damping*sm.fcut_velocity*sample_time);
sm.A_lpf_velocity = 1 + sm.B1_lpf_velocity + sm.B2_lpf_velocity;

%% Noise
sm.sensor_noise_variance = (2e-4)^2*ones(13,1);