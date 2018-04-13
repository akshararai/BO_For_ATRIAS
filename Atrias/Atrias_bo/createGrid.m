function createGrid(grid_slice)
%Slices of size 1000. Change between 1000, 2000, 10000 as starting point by
%changing the grid slice. To change the size of the slice, change n. n=1
%means a 1000 sized slice. To bring slice size down to 100s, change N and
%i_end.
clc; close all; warning off all;
% grid_slice = grid_slice+100;
mkdir(['tmp',num2str(grid_slice)])
cd(['tmp',num2str(grid_slice)])
addpath(genpath('../'))
N = grid_slice*1000;
i_start = N+1;
i_end = N+1000;
disp(['scores_',num2str(i_start),'_',num2str(i_end),'.mat'])

%% model
model_char = 'original';
computer = 'cluster';
model       = strcat('atrias_',model_char);
stopTime = 30.0;
dim = 50;
siminit_name = strcat(['SimInitial/SimInitial_',model_char,'_',computer,'.mat']);
[f1, f2, s, s_4d, simout_vec] = main_ATRIAS_cluster(model, stopTime, dim, siminit_name);
fprintf("costs 1 and 2 - %f and %f, score = %f, score_4d", f1, f2, s); s_4d

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

z = z.xgrid;
z = bsxfun(@times, z, xLimitHigh - xLimitLow);
z = bsxfun(@plus, z, xLimitLow);


%% Evaluate the first 100,000 points
xs = z(i_start:i_end,:);
disp(['points from ',num2str(i_start), ' to ', num2str(i_end)])

%% parpool
% parpool 3
score_grid = zeros(size(xs,1),1);
score_grid_4d = zeros(size(xs,1),4);
fitness = nan(size(xs,1),2);
simout_grid = cell(size(fitness));
sim_t = 2.3410;
stopTime = sim_t + 5;

for iter = 1:size(xs,1)
    [fit1,fit2,score,score_4d,~] = kernel_scores_NM_mismatch(xs(iter,:)', stopTime, model, @findCost, siminit_name);
    fitness(iter,:) = [fit1 fit2];
    score_grid(iter) = score;
    score_grid_4d(iter,:) = score_4d;
    simout_grid{iter} = simout_vec;
%     disp(['points from ',num2str(i_start), ' to ', num2str(i_end)])
    fprintf("iter %f costs 1 and 2 - %f and %f, score = %f, score_4d", iter, fit1, fit2, score); score_4d

end

%%
save(['../data/',model,'scores_4d_',num2str(i_start),'_',num2str(i_end),'_',num2str(dim),'d.mat'], 'score_grid_4d')
save(['../data/',model,'scores_',num2str(i_start),'_',num2str(i_end),'_',num2str(dim),'d.mat'], 'score_grid')
save(['../data/',model,'fitness_',num2str(i_start),'_',num2str(i_end),'_',num2str(dim),'d.mat'], 'fitness')
save(['../data/',model,'simout_',num2str(i_start),'_',num2str(i_end),'_',num2str(dim),'d.mat'], 'simout_grid')
end
