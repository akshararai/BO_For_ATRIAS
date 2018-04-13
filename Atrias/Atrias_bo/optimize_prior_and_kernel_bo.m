function optimize_prior_and_kernel_bo(grun)
grun = grun + 1;
%hwadjust = 0;
hwadjust = 1; % use k^{v_2}_{DoG_{adjust}}
oracle_mean = true;
fake = false;
%% Initialize run-related variables.
warning off all
% Example command to launch a few runs:
% matlab -nodisplay -r "grun=1:5; disturb=1; hwadjust=1; cost_op_str='findCost_metabolic'; optimize_prior_bo"
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
sim_model_char = 'original';
%sim_model_char = 'gear_dynamics';
%sim_model_char = 'without_boom';
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

%% Load simulation-based DoG scores and construct the kernel from them.
load(['data/data_comparisons/atrias_',sim_model_char,'score_1_',num2str(dim),'d.mat']);
dog_scores = scores_1d;
dog_scores(dog_scores < -10) = -10;
if hwadjust >= 1
    opt.hwadjust = hwadjust;
    opt.hwadjust_orig_sim_dog_scores = dog_scores;
    dog_scores_map = [dog_scores, zeros(opt.grid_size,1)];  % sim&hw match initially
end
if hwadjust == 2
    opt.hwadjust_scale = 1; % no scaling for hw dog scores
    locomap(NaN, NaN, opt.grid, dog_scores_map, true);  % init locomap, out_diff=1
    opt.kernel_dims = 1;  % \phi(x)=\phi_{sim}(x)-\bar{g}_{new}
else
    locomap(NaN, NaN, opt.grid, dog_scores_map);  % init locomap
    opt.kernel_dims = size(dog_scores_map,2);
end
ddummy = @(q) 1;  % TODO: handle derivatives better
opt.covfunc = {@covWarp,{@covSEard},@locomap,ddummy,opt.kernel_dims};

%% Load costs collected in simulation and construct a prior for the mean.
% We decided to use costs from 200K short simulations for this version.
num_prior_pts = 200000;
sim_model = strcat('atrias_', sim_model_char);
fitnesses = load(['data/data_comparisons/',sim_model,'fitness_',num2str(dim),'d.mat']);
% 2nd column has findCost_metabolic (1st column has findCost)
sim_costs_all = fitnesses.fitness(:,2);
assert(size(sim_costs_all,1) == num_prior_pts);  % sanity check
if oracle_mean
    mmap(NaN, NaN, opt.grid, sim_costs_all);
else
    gp_prior_mean = load_gp_prior(opt.grid, sim_costs_all, ...
        1:num_prior_pts, num_prior_pts, dim);
    mmap(NaN, NaN, opt.grid, gp_prior_mean);
end
opt.meanfunc = {@meanMmap};  % using mmap() to get sim-based costs

%% BO iterations
num_rand_pts = 2;  % 2 points seems safe for 50D
botrace_prior_and_kernel = {};
for i = 1:length(grun)
    fprintf('------ Prior and Kernel BO run %d ------\n', grun(i));
    rng(grun(i));
    [min_sample, min_value, botrace] = bayesoptBARE(...
        f, opt, model, grun(i), cost_op, num_rand_pts, @locomap, siminit_name);
    botrace_prior_and_kernel{i} = botrace;
end
%go back to initial folder
cd(cwd)
%% Save botrace data
save(char(strcat('data/data_botrace/',sim_model,...
    'botrace_prior_and_kernel_dim_',num2str(dim),...
    '_botrace' , extras_str, '_runs', string(grun(1)), '_', ...
    string(grun(end)),'.mat')), 'botrace_prior_and_kernel','-v7.3')

end

function [gp_prior_mean] = load_gp_prior(grid_all, sim_costs_all, ...
    prior_rand_ids, num_prior_pts, dim)
    % Select num_prior_pts from 50K costs randomly.
    sim_prior_costs = sim_costs_all(prior_rand_ids,:);
    sim_prior_grid = grid_all(prior_rand_ids,:);

    % Simple approach for BO with prior points on opt.grid
    % Not using this, because we need to load 35K pts from 30s simulations,
    % and not all of 200K pts from 5s simulations.
    % locomap(NaN, NaN, opt.grid, sim_costs);
    % opt.meanfunc = {@meanMmap};  % using mmap() to get sim-based costs

    % Construct prior mean from GP. Can work even with off-grid BO.
    % From GPML documentation:
    % A sparse approximate covariance function is defined by composing the
    % apxSparse with a target covariance function and a set of inducing inputs.
    % An inference method is defined by concatenating the struct('s', 0.5)
    % to the infGaussLik inference method: values 0<s<1 -> SPEP
    % (Sparse Power Expectation Propagation).
    fprintf('Loading %d pts into prior mean GP\n', size(sim_prior_costs,1))
    u = sim_prior_grid(1:int32(0.01*num_prior_pts),:);  % inducing points 
    prior_hyp = struct('mean', zeros(1,1), 'cov', zeros(dim+1,1), 'lik', log(0.1));
    prior_hyp.xu = u;  % optimize inducing points
    prior_meanfunc = @meanConst;
    prior_covfunc = {@apxSparse, {@covSEard}, u};
    prior_inffunc = @(varargin) infGaussLik(varargin{:}, struct('s', 0.5));
    prior_likfunc = @likGauss;
    %{
    fprintf('Optimizing prior_hyp...\n')
    prior_hyp2 = minimize(prior_hyp, @gp, -100, prior_inffunc, ...
        prior_meanfunc, prior_covfunc, prior_likfunc, ...
        sim_prior_grid, sim_prior_costs);
    fprintf('Optimized prior_hyp2\n')
    disp(prior_hyp2)
    disp(prior_hyp2.cov)
    %}
    meanfunc = {@meanGP, prior_hyp, prior_inffunc, ...
        prior_meanfunc, prior_covfunc, prior_likfunc, ...
        sim_prior_grid, sim_prior_costs};
    % To speed up BO: evaluate the approximate GP for all grid points, then
    % use the evaluations during BO directly, instead of re-evaluating prior GP.
    gp_prior_mean = feval(meanfunc{:}, prior_hyp, grid_all);
    save(strcat('/tmp/gp_prior_mean.mat'), 'gp_prior_mean');
    fprintf('Done loading sim prior\n')
end
