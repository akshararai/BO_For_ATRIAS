function optimize_kernel_bo(grun)
grun = grun + 1;
hwadjust = 0;
%hwadjust = 1; % use k^{v_2}_{DoG_{adjust}}
fake = 0;
%% Initialize run-related variables.
warning off all
% Example command to launch a few runs:
% matlab -nodisplay -r "grun=1:5; disturb=1; hwadjust=1; cost_op_str='findCost_metabolic'; optimize_kernel_bo"
if ~exist('grun', 'var'), grun=1; end
if ~exist('disturb', 'var'), disturb=0; end
if ~exist('hwadjust', 'var'), hwadjust=0; end
if ~exist('cost_op_str', 'var'), cost_op_str='findCost_metabolic'; end
if ~exist('max_iters', 'var'), max_iters=50; end

extras_str = '';
if disturb, extras_str = char(strcat(extras_str, '_disturb')); end
if hwadjust, extras_str = char(strcat(extras_str, '_hwadjust', num2str(hwadjust))); end

%% Model name and dimensionality

model_char = 'original';
dog_model_char = 'original';
model = strcat('atrias_', model_char);
stopTime = 5.0;
dim = 50;

mkdir(['tmp_kernel_bo',num2str(grun(1)),'_', num2str(grun(end)),extras_str])
cwd = pwd;
cd(['tmp_kernel_bo',num2str(grun(1)),'_', num2str(grun(end)),extras_str])
addpath(genpath('../'))
addpath(genpath(['../../../bayesopt']))
rmpath(genpath('../../../bayesopt/gpml'));  % TODO(rika): remove after update
addpath(genpath('../../../bayesopt/private/gpml_2017'));  % TODO(rika): remove

siminit_name = strcat(['SimInitial_',model_char,'_tolstoy.mat']);
if ~fake
	[f1, f2, s, s_4d, simout_vec] = main_ATRIAS_cluster(model, stopTime, dim, siminit_name);
	fprintf("costs 1 and 2 - %f and %f, score = %f, score_4d", f1, f2, s); s_4d
end

%% Set BO options.

f = str2func('kernel_scores_NM_mismatch');
if fake, f = str2func('fake_kernel_scores_NM_mismatch'); end
cost_op = str2func(cost_op_str);
opt = defaultopt;
opt.optimize_ei = 0;  % TODO: we should use optimize_ei for high-dim BO!
opt.dims = dim;

opt.max_iters = max_iters;
opt.hwadjust_signal_threshold = 50;
opt.hwadjust_scale = 10; % hyperparameter to weigh down the importance of mismatch terms - ideally should be learned
%%
if(dim == 5)
    z = load('grid_5d.mat');
    xLimitHigh = load('xLimitHigh5d.mat');
    xLimitHigh = xLimitHigh.xLimitHigh;
    xLimitLow = load('xLimitLow5d.mat');
    xLimitLow = xLimitLow.xLimitLow;
elseif (dim == 9)
    z = load('grid_9d.mat');
    xLimitHigh = load('xLimitHigh9d.mat');
    xLimitHigh = xLimitHigh.xLimitHigh;
    xLimitLow = load('xLimitLow9d.mat');
    xLimitLow = xLimitLow.xLimitLow;
else 
    z = load('grid_50d.mat');
    xLimitHigh = 1.25*ones(1,50);
    xLimitLow = 0.75*ones(1,50);
end
if opt.optimize_ei == 0
    opt.grid = z.xgrid(1:200000,:);
else
    opt.grid = z.xgrid(1:10000,:);
end
opt.grid_size = size(opt.grid,1);
opt.mins = xLimitLow;
opt.maxes = xLimitHigh;

% Adjust dog_score computed from a hardware-like run to be comparable to
% simulation-based score obatined from 5s runs.
opt.stopTime = 30;
opt.sim_stopTime = 5;  % 10s simulations
num_rand_pts = 2;  % 2 points seems safe for 50D

%% Load simulation-based DoG scores and construct the kernel from them.
tic
locomapF = str2func('locomap');
dog_model = strcat('atrias_', dog_model_char);
load(['data/data_comparisons/',dog_model,'score_1_',num2str(dim),'d.mat']);
dog_scores = scores_1d;
dog_scores(dog_scores < -10) = -10;
if hwadjust >= 1
    opt.hwadjust = hwadjust;
    opt.hwadjust_orig_sim_dog_scores = dog_scores;
    dog_scores_map = [dog_scores, zeros(opt.grid_size,1)];  % sim&hw match initially
else
    dog_scores_map = dog_scores;
end
% scoresMap = containers.Map;
% scoresDim = size(dog_scores,2);
% for i = 1:opt.grid_size
%     mkey = sprintf('%0.4g', opt.grid(i,:));
%     scoresMap(mkey) = dog_scores(i,:);
% end
if hwadjust == 2
    opt.hwadjust_scale = 1; % no scaling for hw dog scores
    locomapF(NaN, NaN, opt.grid, dog_scores_map, true);  % init locomap, out_diff=1
    opt.kernel_dims = 1;  % \phi(x)=\phi_{sim}(x)-\bar{g}_{new}
else
    locomapF(NaN, NaN, opt.grid, dog_scores_map);  % init locomap
    opt.kernel_dims = size(dog_scores_map,2);
end
ddummy = @(q) 1;  % TODO: handle derivatives better
opt.covfunc = {@covWarp,{@covSEard},locomapF,ddummy,opt.kernel_dims};

%% BO iterations
botrace_kernel = {};
for i = 1:length(grun)
    fprintf('------ DoG kernel BO run %d ------\n', grun(i));
    rng(grun(i));
    [min_sample, min_value, botrace] = bayesoptBARE(...
        f, opt, model, cost_op, num_rand_pts, locomapF, siminit_name);
    botrace_kernel{i} = botrace;
end
toc
%go back to initial folder
cd(cwd)
%% Save botrace data
save(char(strcat('data/data_botrace/',dog_model,'botrace_kernel_dim_',num2str(dim),'_botrace' , extras_str, ...
    '_runs', string(grun(1)), '_', string(grun(end)),'.mat')), ...
    'botrace_kernel','-v7.3')

end
