%% Avoid virtual leg length limits by ensuring we don't command position/velocities/accelerations that would hit them.
function [l_leg_cmd, dl_leg_cmd, ddl_leg_cmd] ...
    = AvoidVirtualLengthLimits(l_leg_desired, dl_leg_desired, ddl_leg_desired, ...
                               min_leg_length, max_leg_length)

    min_length_distance = l_leg_desired - min_leg_length;
    max_length_distance = max_leg_length - l_leg_desired;

    l_leg_cmd = min(max(l_leg_desired, min_leg_length), max_leg_length);
    dl_leg_cmd = dl_leg_desired;
    ddl_leg_cmd = ddl_leg_desired;

    distance_threshold = 0.02;

    if (min_length_distance < distance_threshold && sign(dl_leg_desired) < 0) || ...
       (max_length_distance < distance_threshold && sign(dl_leg_desired) > 0)
        dl_leg_cmd = 0;
    end
    
    if (min_length_distance < distance_threshold && sign(ddl_leg_desired) < 0) || ...
       (max_length_distance < distance_threshold && sign(ddl_leg_desired) > 0)
        ddl_leg_cmd = 0;
    end

end