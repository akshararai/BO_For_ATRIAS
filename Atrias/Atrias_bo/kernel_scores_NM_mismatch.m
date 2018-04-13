function [fit1, fit2, score, score_4d, simout_vec] = kernel_scores_NM_mismatch(param, stopTime, model, cost_op, siminitial_name)
if(nargin < 4)
    cost_op = [];
end
if(nargin < 5)
    siminitial_name = [];
end

dim = length(param);
%%
    
    nan_val = cost_op();
    score = 0;
    fit1 = nan_val;
    fit2 = nan_val;
    score_4d = [0 0 0 0];
    simout_vec = [];
    
%% Assign
    if(dim ~= 50)
        nmm_param = ones(50,1);
    else
        nmm_param = param;
    end
    [sm, update_freq, sample_time] = fn_set_SimParams();
    sm = fn_set_MechParams(sm);    nm  = fn_set_vnmc_params(sm, sample_time);
    [nm, nmc] = fn_set_vnmc_CtrlParams(nm, nmm_param);
    nmc_BLCtrlPar   = nmc.BLCtrlPar;
    nmc_FPCtrlPar   = nmc.FPCtrlPar;
    nmc_LLCtrlPar   = nmc.LLCtrlPar;
    nmc_StCtrlPar   = nmc.StCtrlPar;
    nmc_SwCtrlPar   = nmc.SwCtrlPar;
    nmc_TrCtrlPar   = nmc.TrCtrlPar;

    nm1.l_v     = nm.l_v;
    nm1.FmaxVAS = nm.FmaxVAS;
    assignin('base', 'nm', nm);
    assignin('base', 'nm1', nm1);
    assignin('base', 'nmc_BLCtrlPar', nmc_BLCtrlPar);
    assignin('base', 'nmc_FPCtrlPar', nmc_FPCtrlPar);
    assignin('base', 'nmc_LLCtrlPar', nmc_LLCtrlPar);
    assignin('base', 'nmc_StCtrlPar', nmc_StCtrlPar);
    assignin('base', 'nmc_SwCtrlPar', nmc_SwCtrlPar);
    assignin('base', 'nmc_TrCtrlPar', nmc_TrCtrlPar);

    [mc, idc, rc] = set_control_params_raibert_id(sm, sample_time);
    if(dim==5)
        rc_CtrlPar = [param, rc.desired_theta, rc.desired_z, rc.kpz, rc.kdz];
    elseif(dim==9)
        rc_CtrlPar = param;
    else
        rc_CtrlPar = [rc.k_placement, rc.T_step, rc.Cd, rc.kpt, rc.kdt, rc.desired_theta, rc.desired_z, rc.kpz, rc.kdz];
    end
    assignin('base','rc_CtrlPar', rc_CtrlPar);

    
    %% Simulate
    if(dim==50)
    load(siminitial_name);
    assignin('base', 'SimInitial', SimInitial); 
    set_param(model,'LoadInitialState','on','InitialState','SimInitial');
    end

    try
    simout = sim(model,...
                'SimulationMode','accelerator',...
                'StopTime', num2str(stopTime),  ...
                'SaveOutput',                       'on');
    catch exception
        exception.message
        f = nan;
        return;
    end
    set_param(model,'LoadInitialState','off');

    %% Evaluate score

    %simulate each sample and store cost
    %parfor i = 1:popSize
%     for i = 1:popSize
    [fit1, fit2, score, score_4d] = evaluateGaitScore(simout, stopTime, cost_op);
    simout_vec = simout_to_vec(simout);%     end
    if(isnan(score))
        score = 0;
    end
    if(isnan(fit1))
        fit1 = nan_val;
    end
    if(isnan(fit2))
        fit2 = nan_val;
    end
    