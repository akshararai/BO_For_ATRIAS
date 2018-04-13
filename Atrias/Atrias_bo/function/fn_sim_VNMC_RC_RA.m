function f = fn_sim_VNMC_RC_RA(model, rtp, nm, param, v_tgt, c_weight, flag_disp)

c_d     = c_weight{1};
c_v     = c_weight{2};
c_tau	= c_weight{3};

[nm, nmc, parOffR] = fn_set_vnmc_CtrlParams(nm, param);
nmc_BLCtrlPar   = nmc.BLCtrlPar;
nmc_FPCtrlPar   = nmc.FPCtrlPar;
nmc_LLCtrlPar   = nmc.LLCtrlPar;
nmc_StCtrlPar   = nmc.StCtrlPar;
nmc_SwCtrlPar   = nmc.SwCtrlPar;
nmc_TrCtrlPar   = nmc.TrCtrlPar;

paramSets = Simulink.BlockDiagram.modifyTunableParameters(rtp, ...
            'nmc_BLCtrlPar',            nmc_BLCtrlPar, ...
            'nmc_FPCtrlPar',            nmc_FPCtrlPar, ...
            'nmc_LLCtrlPar',            nmc_LLCtrlPar, ...
            'nmc_StCtrlPar',            nmc_StCtrlPar, ...
            'nmc_SwCtrlPar',            nmc_SwCtrlPar, ...
            'nmc_TrCtrlPar',            nmc_TrCtrlPar);

if flag_disp
    tic;
end
try
    simout = sim(model,...
            'SimulationMode','rapid',...
            'RapidAcceleratorParameterSets',    paramSets, ...
            'RapidAcceleratorUpToDateCheck',    'off', ...
            'SaveOutput',                       'on');
catch exception
%     exception.message
    f = nan;
    return;
end
sim_tout    = get(simout, 'tout');
sim_tEnd    = sim_tout(end);
if flag_disp
    tend    = toc;
    treal = sim_tEnd/tend*100;
    disp(['simulation speed: ' num2str(treal) '% of real time']);
end
        
simout_flag             = get(simout, 'simout_flag');
simout_xTorso_switch	= get(simout, 'simout_xTorso_switch');
% simout_yTorso           = get(simout, 'simout_yTorso');
% simout_nStep            = get(simout, 'simout_nStep');

t_switch	= simout_xTorso_switch.time;

flag = simout_flag.signals.values;
if flag == 1
    simout_Tau2         = get(simout, 'simout_Tau2');
    
    tau2    = simout_Tau2.signals.values;
    dt      = sim_tEnd - t_switch;
    
    f = c_tau*tau2/dt + c_v*(v_ave - v_tgt)^2;
else
    simout_xTorsoEnd	= get(simout, 'simout_xTorsoEnd');
    simout_dxTorso      = get(simout, 'simout_dxTorso');
    
    d       = simout_xTorsoEnd.signals.values;
        ii  = find(simout_dxTorso.time > t_switch, 1);
    v_ave	= mean(simout_dxTorso.signals.values(ii:end));
    
    f = 1e6 - c_d*d + c_v*(v_ave - v_tgt)^2;
end

f = f + parOffR;
