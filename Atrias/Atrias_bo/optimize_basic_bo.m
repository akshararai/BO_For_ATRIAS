function optimize_basic_bo(grun)
grun = grun + 1; %grun goes from 0 to N-1. Make MATLAB compatible
fake = false; % turn to true if you want to manually enter costs, avoid sim
%% Initialize run-related variables.
warning off all
% Example command to launch a few runs:
% matlab -nodisplay -r "grun=1:5; disturb=1; cost_op_str='findCost_metabolic'; optimize_prior_bo"
%Defaults:
if ~exist('grun', 'var'), grun=1; end 
if ~exist('cost_op_str', 'var'), cost_op_str='findCost_metabolic'; end
if ~exist('max_iters', 'var'), max_iters=50; end

extras_str = '';

%% Model name and dimensionality

model_char = 'original'; %SE is always run on the original
model = strcat('atrias_', model_char);
dim = 5;

%create a new directory for each grun and move into it. 
%This is to avoid simulink errors on running in parallel on the cluster
mkdir(['tmp_basic_bo',num2str(grun(1)),'_', num2str(grun(end)),extras_str])
cwd = pwd;
cd(['tmp_basic_bo',num2str(grun(1)),'_', num2str(grun(end)),extras_str])
addpath(genpath('../'))
addpath(genpath(['../../../bayesopt']))

%siminit is needed in 50d. Ask Akshara on help to create and use it.
%TODO Akshara: Change code to not need siminit for simple version.
siminit_name = strcat(['SimInitial_',model_char,'_cluster.mat']);

if ~fake %initialize atrias model variables, run a short simulation to make sure everything is initialized, etc.
    tmp_stopTime  = 5.0;
    [f1, f2, s, s_4d, simout_vec] = main_ATRIAS_cluster(model, tmp_stopTime, dim, siminit_name);
    fprintf("costs 1 and 2 - %f and %f, score = %f, score_4d", f1, f2, s); s_4d
end

%% Set BO options. 

f = str2func('kernel_scores_NM_mismatch');
if fake, f = str2func('fake_kernel_scores_NM_mismatch'); end
cost_op = str2func(cost_op_str);
opt = defaultopt;
opt.optimize_ei = 0;  % TODO: we should use optimize_ei for high-dim BO!
opt.dims = dim;
opt.covfunc = {@covSEard};
opt.max_iters = max_iters;
%% Load the search limits for the controller space
if(dim == 5)
    z = load('grid_5d.mat');
    xLimitHigh = load('xLimitHigh5d.mat');
    xLimitHigh = xLimitHigh.xLimitHigh;
    xLimitLow = load('xLimitLow5d.mat');
    xLimitLow = xLimitLow.xLimitLow;
    num_grid_points = 20000;
elseif (dim == 9)
    z = load('grid_9d.mat');
    xLimitHigh = load('xLimitHigh9d.mat');
    xLimitHigh = xLimitHigh.xLimitHigh;
    xLimitLow = load('xLimitLow9d.mat');
    xLimitLow = xLimitLow.xLimitLow;
    num_grid_points = 100000;
else 
    z = load('grid_50d.mat');
    xLimitHigh = 1.25*ones(1,50);
    xLimitLow = 0.75*ones(1,50);
    num_grid_points = 200000;
end
opt.grid = z.xgrid(1:num_grid_points,:);
opt.grid_size = size(opt.grid,1);
opt.mins = xLimitLow;
opt.maxes = xLimitHigh;

% Adjust dog_score computed from a hardware-like run to be comparable to
% simulation-based score obatined from 5s runs.
opt.stopTime = 30;
opt.sim_stopTime = 5;  % 10s simulations

%% BO iterations
num_rand_pts = 2;  % number of points to choose randomly before BO.
botrace_prior = {};
for i = 1:length(grun)
    fprintf('------ Prior BO run %d ------\n', grun(i));
    rng(grun(i));
    [min_sample, min_value, botrace] = bayesoptBARE(...
        f, opt, model, cost_op, num_rand_pts, [], siminit_name);
    botrace_prior{i} = botrace;
end
%go back to initial folder
cd(cwd)
%% Save botrace data
save(char(strcat('data/data_botrace/',model,...
    'botrace_SEard_dim_',num2str(dim),...
    '_botrace' , extras_str, '_runs', string(grun(1)), '_', ...
    string(grun(end)),'.mat')), 'botrace_prior','-v7.3')

end