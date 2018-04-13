% Create SimInitial.m (used for first steps before NM controller takes over).
warning off all
mkdir('tmp_create_sim_init')
cwd = pwd;

model_char = 'original';
model = strcat('atrias_', model_char);
open(model)
cd('tmp_create_sim_init')
addpath(genpath('../'))
% First stage of generating SimInitial: initialize all values in workspace
[fit1, fit2, score, simout] = main_ATRIAS_cluster(model,5,5)

nStep_thr = 10;
set_param(model,'SaveFinalState','on','FinalStateName',...
'xFinal','SaveCompleteFinalSimState','on');
try
simout = sim(model,...
            'SimulationMode','accelerator',...
            'StopTime', num2str(10),  ...
            'SaveOutput',                       'on');
catch exception
    exception.message
    f = nan;
    return;
end

nStep = get(simout,'simout_nStep');
check_nstep = nStep.signals.values(end);
assert(check_nstep == 10);
sim_t = nStep.time(end) % e.g. 3.7160 for undisturbed SimInitial
nStep_thr = 1000;
test_VNMC = 1;
% Second stage: generate and save SimInitial
set_param(model,'SaveFinalState','on','FinalStateName',...
'xFinal','SaveCompleteFinalSimState','on');
try
simout = sim(model,...
            'SimulationMode','accelerator',...
            'StopTime', num2str(sim_t),  ...
            'SaveOutput',                       'on');
catch exception
    exception.message
    f = nan;
    return;
end

nStep = get(simout,'simout_nStep');
check_nstep = nStep.signals.values(end);
assert(check_nstep == 10);
xFinal=get(simout,'xFinal')
SimInitial=xFinal;
cd(cwd)  % go back to parent folder

siminit_name = strcat(['SimInitial/SimInitial_',model_char,'_tolstoy.mat']);
save(siminit_name,'SimInitial')

stopTime  = sim_t+5;
set_param(model,'LoadInitialState','on','InitialState', 'SimInitial');
try
simout = sim(model,...
            'SimulationMode','accelerator',...
            'StopTime', num2str(stopTime),  ...
            'SaveOutput',                       'on');
catch exception
    exception.message
    return;
end

[fit1, fit2, score, score_4d, simout_vec] = evaluateGaitScore(simout, stopTime)

%{
%% Check a good point.
pt2_unscaled = [...
    0.9309    0.9996    1.2497    1.2032    0.9929    0.9210    1.1385 ...
    0.7793    1.1211    0.7571    0.9437    1.1732    0.8495    1.0013 ...
    1.2395    0.9331    0.9547    1.1367    1.1546    1.2315    0.9703 ...
    1.1634    0.9122    1.1261    1.1247    1.2149    0.7907    1.1290 ...
    1.2180    1.0918    0.8577    1.0266    0.9809    0.9394    1.1416 ...
    1.1253    0.8837    1.1451    0.9774    1.1587    1.1242    0.9356 ...
    0.9080    0.7957    1.0765    0.9358    0.7551    0.9057    1.0435 ...
    0.8036 ];

[cost, ~, dog_score, ~, simout_vec] = kernel_scores_NM(...
    pt2_unscaled, 30, model, str2func('findCost_metabolic'))

assert(cost > 2.5 && cost < 3.5)
assert(dog_score > 120 && dog_score < 130)

%% Check a bad point.
[f1, f2, s, model, simout, stopTime] = main_ATRIAS_cluster(0,0,disturb,NaN);
pt1_unscaled = ones(1,50)*0.5;

[cost, ~, dog_score, ~, simout_vec] = kernel_scores_NM(...
    pt1_unscaled, 30, model, str2func('findCost_metabolic'))

assert(cost > 97 && cost < 99)
assert(dog_score < 1)
%}
