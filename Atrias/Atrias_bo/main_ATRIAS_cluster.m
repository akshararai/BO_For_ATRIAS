function [fit1, fit2, score, score_4d, simout_vec, siminitial_name] = main_ATRIAS_cluster(model, sim_t, dim, siminitial_name,perturb)
%%
fit1 = nan;
fit2 = nan;
score = 0;
score_4d = zeros(4,1);
simout_vec = [];
if(nargin < 5)
    perturb = 0;
end
% siminitial_name = 'SimInitial/SimInitial_original_cluster.mat'; % replace by machine name
%%

addpath(genpath('../function/'));
addpath('../param/')
addpath('../../common/')
addpath('../../images/')
addpath('../../models/')
addpath(genpath('../../bayesopt'));

% s_VNMC_CtrlPar	 = 'param_handtuned';
% load(s_VNMC_CtrlPar)

nStep_switchCtrl        = 10;
test_VNMC = 0;
if(dim==50), test_VNMC = 1; end
nStep_thr = 1000;

assignin('base','test_VNMC', test_VNMC);
assignin('base','nStep_switchCtrl', nStep_switchCtrl);
assignin('base','nStep_thr', nStep_thr);
if dim==50
    load(siminitial_name);
    assignin('base', 'SimInitial', SimInitial); 
end
assignin('base', 'sim_t', sim_t);

[sm, update_freq, sample_time] = fn_set_SimParams();
sm = fn_set_MechParams(sm,perturb);
colors;
%% set VNMC parameters
param = ones(50,1);
nm  = fn_set_vnmc_params(sm, sample_time);
[nm, nmc] = fn_set_vnmc_CtrlParams(nm, param);
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


%% set raibert control parameters

[mc, idc, rc] = set_control_params_raibert_id(sm, sample_time);

assignin('base', 'rc', rc);
assignin('base', 'idc', idc);
assignin('base', 'mc', mc);
assignin('base', 'sm', sm);
assignin('base', 'sample_time', sample_time);

rc_CtrlPar = [rc.k_placement, rc.T_step, rc.Cd, rc.kpt, rc.kdt, rc.desired_theta, rc.desired_z, rc.kpz, rc.kdz];
assignin('base','rc_CtrlPar', rc_CtrlPar);

v_tgt_step= rc.desired_xd*ones(10000,1);
assignin('base', 'v_tgt_step', v_tgt_step)
%% load and simulate
load_system(model);
stopTime = sim_t;
 
if(dim==50)
set_param(model,'LoadInitialState','on','InitialState', 'SimInitial');
end
 
try
simout = sim(model,...
             'SimulationMode','accelerator',...
             'StopTime', num2str(stopTime),  ...
             'SaveOutput',                       'on');
catch exception
     exception.message
     return;
end
%% score
[fit1, fit2, score, score_4d, simout_vec] = evaluateGaitScore(simout, stopTime);
end

