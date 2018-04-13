function [cost] = findCost_metabolic(sim_vars, stopTime, get_op)
    WORST_METABOLIC_COST = 120;
    cost = WORST_METABOLIC_COST;
    if nargin < 1
        return;  % return worst cost
    end
    if nargin < 2
        stopTime = 30;
    end
    if(nargin < 3)
        get_op = str2func('get');  % get for simulink simout by default
    end
    sv = false;

%     VNMC_on = get(sim_vars, 'VNMC_on');
%     vnmc_tmp = VNMC_on.signals.values;
%     vnmc(:) = vnmc_tmp(1,1,:);
%     idxs = find(vnmc);
    start_idx = 1;%idxs(1);
%     diff = idxs(end) - idxs(1);
    time = get_op(sim_vars,'time');
    time = time.signals.values;
   
    if sv
        tmstamp = strcat(datestr(clock,'yyyy-mm-dd-HHMM'),'m',...
            datestr(clock,'ss'),'s');
        time_end = time(end);
        fname = sprintf('/tmp/findCost_mtb_%s_t_%0.4f.mat', tmstamp, time_end);
    end
 
    xTorsoEnd = get_op(sim_vars, 'xTorsoEnd');
    xTorsoEnd = xTorsoEnd.signals.values;
    simout_xTorsoEnd = xTorsoEnd(end);
    
    
    if time(end) < stopTime
        cost = 100 - simout_xTorsoEnd;
        % cost = 50 - time(end);
        if sv, save(fname, 'simout_xTorsoEnd','time_end','cost','-v7.3'); end
        return
    end
    
    LStance = get_op(sim_vars,'LStance');
    LStance = LStance.signals.values;
    RStance = get_op(sim_vars,'RStance');
    RStance = RStance.signals.values;
    
    rsteps = find(RStance);
    lsteps = find(LStance);

    if (length(LStance) < 30 || length(RStance) < 30)
        return;
    end
    if(length(rsteps) < 1 || length(lsteps) < 1)
        return;
    end

    step_idx = [];
    step_idx(1) = min(rsteps(1), lsteps(1));
    doubleStance = LStance & RStance;
    steps = find(doubleStance);
    step_idx(2) = steps(1);
    for i=1:length(steps)-1
        if(steps(i+1)-steps(i)>1)
            step_idx = [step_idx, steps(i+1)];
        end 
    end
    stanceFlag = zeros(size(LStance));
    stanceFlag(step_idx) = 1;
    nSteps = sum(stanceFlag);
    
    % idx2 = floor(find(simout == stopTime)/4);
    
    CoMdx_sim = get_op(sim_vars, 'CoMVelocity');
    CoMdx_sim = CoMdx_sim.signals.values;
    CoMdx_sim = CoMdx_sim(start_idx:end,:);
    
    motor_torques = get_op(sim_vars, 'motor_torques');
    motor_torques = motor_torques.signals.values;
    motor_torques = motor_torques(start_idx:end,:);
    
    motor_cost = sum(sum(motor_torques(10:end,:).^2))*10^-9;

    CoMdx = CoMdx_sim;
    v_tgt = 1.0*ones(size(CoMdx,1), size(CoMdx,2));

    v_avg = zeros(nSteps-1,1);
    v_tgt_step = v_avg;
    for step = 1:nSteps-1
        v_avg(step) = mean(CoMdx(step_idx(step):step_idx(step+1)));
        v_tgt_step(step) = v_tgt(step_idx(step)+1);
    end
    err = sum((v_tgt_step - v_avg).^2);
    cost = err + motor_cost;

    if sv
        save(fname, 'v_tgt','v_tgt_step','v_avg','simout_xTorsoEnd',...
            'time_end','CoMdx','err','motor_torques','motor_cost','cost','-v7.3')
    end
    
end
