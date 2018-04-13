% get_RA(model, sim_t, sim_d, flag_vnmc_rc_hybrid)

colors;     % sm

param   = zeros(50,1);

[sm, sample_freq, sample_time] = fn_set_SimParams(sm);
sm = fn_set_MechParams(sm);

% set VNMC parameters
nm  = fn_set_vnmc_params(sm, sample_time);
[nm, nmc] = fn_set_vnmc_CtrlParams(nm, param);

% set raibert control parameters
[mc, idc, rc] = set_control_params_raibert_id(sm, sample_time);
nmc_BLCtrlPar   = nmc.BLCtrlPar;
nmc_FPCtrlPar   = nmc.FPCtrlPar;
nmc_LLCtrlPar   = nmc.LLCtrlPar;
nmc_StCtrlPar   = nmc.StCtrlPar;
nmc_SwCtrlPar   = nmc.SwCtrlPar;
nmc_TrCtrlPar   = nmc.TrCtrlPar;

load_system(model);
rtp = Simulink.BlockDiagram.buildRapidAcceleratorTarget(model);

nm1.l_v     = nm.l_v;
nm1.FmaxVAS = nm.FmaxVAS;