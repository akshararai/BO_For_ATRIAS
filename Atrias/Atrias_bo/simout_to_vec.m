% Extract information from simuout that is used for cost computation 
% and record these in a vector to store and use for NN training.
function [simout_vec] = simout_to_vec(simout)
    try
        simout_vec = simout_to_vec_raw(simout);
    catch
        disp('simout failed')
        keyboard
        simout_vec = [];  % findCost*() returns worst cost on this
    end
end

function [simout_vec] = simout_to_vec_raw(sim_vars)
    simout_vec = [];
    
    VNMC_on = get(sim_vars, 'VNMC_on');
    vnmc_tmp = VNMC_on.signals.values;
    vnmc(:) = vnmc_tmp(1,1,:);
    vnmc = vnmc';
    simout_vec.vnmc = vnmc(1:10:end,:);
    
    % Fields used by findCost.m    
    time = get(sim_vars,'time');
    time = time.signals.values;
%     simout_vec(1) = time(end);
%     simout_vec(2) = length(time);
    simout_vec.time = time(1:10:end,:);
    
    xTorsoEnd = get(sim_vars, 'xTorsoEnd');
    xTorsoEnd = xTorsoEnd.signals.values;
%     simout_xTorsoEnd = xTorsoEnd(end);
%     simout_vec(3) = simout_xTorsoEnd;
    simout_vec.xTorso = xTorsoEnd(1:10:end,:);
    
    CoMdx_sim = get(sim_vars, 'CoMVelocity');
    CoMdx = CoMdx_sim.signals.values;
%     errVelocity = sum((dx_des - CoMdx).^2);
%     simout_vec(4) = errVelocity;
    simout_vec.CoMdx = CoMdx(1:10:end,:);
    
    simout_Theta = get(sim_vars, 'TorsoPitch');
    theta_tmp = simout_Theta.signals.values;
    theta(:) = theta_tmp(:);
    theta = theta';
%     simout_vec(6) = theta(end);
    simout_vec.theta = theta(1:10:end,:);
    
%     Velocity_cost = get(sim_vars, 'Velocity_cost');
%     vel_tmp = Velocity_cost.signals.values;
%     vel = vel_tmp(:);
%     simout_vec(7) = vel(end);

    simout_Tau = get(sim_vars, 'motor_torques');
    tau = simout_Tau.signals.values;
%     simout_vec(8) = tau;
    simout_vec.tau = tau(1:10:end,:);
    
    swingFootz = get(sim_vars, 'PrimaryToSecondaryZ');
    swingFootz = swingFootz.signals.values;
    simout_vec.swingFootz = swingFootz(1:10:end,:);
    
    LStance = get(sim_vars, 'LStance');
    LStance = LStance.signals.values;
    simout_vec.LStance = LStance(1:10:end,:);
    
    RStance = get(sim_vars, 'RStance');
    RStance = RStance.signals.values;
    simout_vec.RStance = RStance(1:10:end,:);
    
    comz = get(sim_vars, 'CoMHeight');
    comz = comz.signals.values;
    simout_vec.comz = comz(1:10:end,:);
    
    EnteringStance = get(sim_vars, 'EnteringStance');
    EnteringStance = EnteringStance.signals.values;
    simout_vec.EnteringStance = EnteringStance(1:10:end,:);
end
