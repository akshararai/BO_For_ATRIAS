function optimize_nn_bo(grun)
grun = grun + 1;
%nn_dim = 4;  % 4D or 13D simout vecs reconstructed by NN
nn_dim = 13;  % 4D or 13D simout vecs reconstructed by NN
fake = false;
cost_op_str='findCost_smooth';

%% Initialize run-related variables.
warning off all
% Example command to launch a few runs:
% matlab -nodisplay -r "grun=1:5; disturb=1; cost_op_str='findCost_metabolic'; optimize_kernel_bo"
if ~exist('grun', 'var'), grun=1; end
if ~exist('disturb', 'var'), disturb=0; end
if ~exist('cost_op_str', 'var'), cost_op_str='findCost_metabolic'; end
if ~exist('max_iters', 'var'), max_iters=50; end

extras_str = '';
if disturb, extras_str = char(strcat(extras_str, '_disturb')); end

%% Model name and dimensionality

model_char = 'original';
sim_model_char = 'original';
model = strcat('atrias_', model_char);
stopTime = 5.0;
dim = 50;

mkdir(['tmp_nn_bo',num2str(grun(1)),'_', num2str(grun(end)),extras_str])
cwd = pwd;
cd(['tmp_nn_bo',num2str(grun(1)),'_', num2str(grun(end)),extras_str])
addpath(genpath('../'))
addpath(genpath(['../../../bayesopt']))
rmpath(genpath('../../../bayesopt/gpml'));  % TODO(rika): remove after update
addpath(genpath('../../../bayesopt/private/gpml_2017'));  % TODO(rika): remove

siminit_name = strcat(['SimInitial_',model_char,'_cluster.mat']);

if ~fake
    [f1, f2, s, s_4d, simout_vec] = main_ATRIAS_cluster(model, stopTime, dim, siminit_name);
    fprintf("costs 1 and 2 - %f and %f, score = %f, score_4d", f1, f2, s); s_4d
end


%% Set BO options.

f = str2func('kernel_scores_NM_mismatch');
if fake, f = str2func('fake_kernel_scores_NM_mismatch'); end
cost_op = str2func(cost_op_str);
opt = defaultopt;
opt.dims = dim;
opt.kernel_dims = nn_dim;

opt.max_iters = max_iters;
% opt.hwadjust_signal_threshold = 50;
%%
if(dim == 5)
    z = load('grid_5d.mat');
    xLimitHigh = load('xLimitHigh5d.mat');
    xLimitHigh = xLimitHigh.xLimitHigh;
    xLimitLow = load('xLimitLow5d.mat');
    xLimitLow = xLimitLow.xLimitLow;
    opt.grid = z.xgrid(1:end,:);
elseif (dim == 9)
    z = load('grid_9d.mat');
    xLimitHigh = load('xLimitHigh9d.mat');
    xLimitHigh = xLimitHigh.xLimitHigh;
    xLimitLow = load('xLimitLow9d.mat');
    xLimitLow = xLimitLow.xLimitLow;
    opt.grid = z.xgrid(1:end,:);
else 
    z = load('grid_50d.mat');
    xLimitHigh = 1.25*ones(1,50);
    xLimitLow = 0.75*ones(1,50);
    opt.grid = z.xgrid(1:200000,:);
end
opt.grid_size = size(opt.grid,1);
opt.mins = xLimitLow;
opt.maxes = xLimitHigh;

% Adjust dog_score computed from a hardware-like run to be comparable to
% simulation-based score obatined from 5s runs.
opt.stopTime = 30;
opt.sim_stopTime = 5;  % 10s simulations

%% Create DoG scores map.
% TODO(akshara): is the data in parent dir on the cluster or not?
nn_fname = strcat('../data/data_comparisons/','nn_out_vals_',sim_model_char,...
    '_200K_from',num2str(nn_dim),'Dsimouts.mat');
fprintf('loading NN from %s\n', nn_fname);
assert(exist(nn_fname, 'file') == 2)
load(nn_fname);  % nn_out_vals
assert(nn_dim == size(nn_out_vals,2));
nn_out_vals = nn_out_vals(1:opt.grid_size,:);

locomap(NaN, NaN, opt.grid, nn_out_vals);  % init locomap
ddummy = @(q) 1;  % TODO: handle derivatives better if possible
opt.covfunc = {@covWarp,{@covSEard},@locomap,ddummy,nn_dim};

%% BO iterations
num_rand_pts = 2;  % 2 points seems safe for 50D
botrace_nn = {};
for i = 1:length(grun)
    fprintf('------ NN%d kernel BO in %dd run %d ------\n', nn_dim, dim, grun(i));
    rng(grun(i));
    [min_sample, min_value, botrace] = bayesoptBARE(...
        f, opt, model, grun(i), cost_op, num_rand_pts, @locomap, siminit_name);
    botrace_nn{i} = botrace;
end
%go back to initial folder
cd(cwd)
%% Save botrace data
save(char(strcat('data/data_botrace/atrias_',sim_model_char,...
    'botrace_nn',num2str(nn_dim),'_dim_',num2str(dim),'_botrace' , extras_str, ...
    '_runs', string(grun(1)), '_', string(grun(end)),'.mat')), ...
    'botrace_nn','-v7.3')

end
