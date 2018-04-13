function [l_leg, dl_leg, ddl_leg] = LegLengthTrajectory(t, l_leg_initial, l_leg_final, t_trajectory, max_velocity, max_acceleration)
    t = max(0,t);    
    t_trajectory = max(0, t_trajectory);
    target_distance = l_leg_final - l_leg_initial;
    [t_ramp, t_hold] = TrapezoidalTrajectoryTimings(t_trajectory, max_acceleration, max_velocity, abs(target_distance));
    max_acceleration = max_acceleration*sign(target_distance);

    % generate trajectory
    if (t >= 0 && t <= t_ramp)
        l_leg = l_leg_initial + 0.5*max_acceleration*t^2;
        dl_leg = max_acceleration*t;
        ddl_leg = max_acceleration;
    elseif (t > t_ramp && t <= t_ramp+t_hold)
        l_leg = l_leg_initial + 0.5*max_acceleration*t_ramp^2 + (t-t_ramp)*max_acceleration*t_ramp;
        dl_leg = max_acceleration*t_ramp;
        ddl_leg = 0;
    elseif (t > t_ramp+t_hold && t <= t_hold + 2*t_ramp)
        l_leg = l_leg_initial + 0.5*max_acceleration*t_ramp^2 + t_hold*max_acceleration*t_ramp + ...
                max_acceleration*t_ramp*(t-t_ramp-t_hold) - 0.5*max_acceleration*(t-t_ramp-t_hold)^2;
        dl_leg = max_acceleration*t_ramp - max_acceleration*(t-t_ramp-t_hold);
        ddl_leg = -max_acceleration;
    else
        l_leg = l_leg_initial + max_acceleration*t_ramp^2 + max_acceleration*t_ramp*t_hold;
        dl_leg = 0;
        ddl_leg = 0;
    end
end

