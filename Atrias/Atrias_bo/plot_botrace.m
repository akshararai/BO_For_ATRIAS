
%% Analyze explicit mismatch simulation experiments.

% Plots for BO results with error region.
addpath(genpath('../../tools'))  % for plot_botrace_helper

%% Load BO runs data.
dir_prefix = '../../../boforloco_runs_writeups/runs_mismatch/';
cost_max = 1.35; cost_good = 0.7; failed_cost = 2; cost_min = 0.65;
%cost_max = 100; cost_good = 20; failed_cost = 120; cost_min = 0;

%{
algodirs={...
    'SE',...
    'nn',...
    'kernel',...
    'hwadjust',...
    };
algofiles={...
    'originalbotrace_SEard_50',...
    'gear_dynamicsbotrace_nn13_dim_50',...
    'gear_dynamicsbotrace_kernel_50',...
    'gear_dynamicsbotrace_kernel_dim_50hwadjust1_',...
    };
clrs={'r', 'g', 'k', 'b'};
markers={'none', '^', '*', '+'};
lstyles={'-', '-', '-', '-'};
labels={'SE kernel', 'trajNN knl', 'DoG knl', 'adj DoG knl'};
%}
%{
algodirs={...
    'SE',...
    'nn',...
    'kernel',...
    'hwadjust'
    };
algofiles={...
    'originalbotrace_SEard_50',...
    'without_boombotrace_nn13_dim_50',...
    'without_boombotrace_kernel_50',...
    'without_boombotrace_kernel_dim_50hwadjust1_'
    };
clrs={'r', 'g', 'k', 'b'};
markers={'none', '^', '*', '+'};
lstyles={'-', '-', '-', '-'};
labels={'SE kernel', 'trajNN knl', 'DoG knl', 'adj DoG knl'};
%}

% Smooth cost plots.
%{
algodirs={...
    'SE',...
    'nn',...
    'kernel',...
    'hwadjust'
    };
algofiles={...
    'originalbotrace_SEard_dim_50smooth',...
    'gear_dynamicsbotrace_nn13_dim_50smooth',...
    'gear_dynamicsbotrace_kernel_dim_50smooth_',...
    'gear_dynamicsbotrace_kernel_dim_50smooth_hwadjust1_'
    };
clrs={'r', 'g', 'k', 'b'};
markers={'none', '^', '*', '+'};
lstyles={'-', '-', '-', '-'};
labels={'SE kernel', 'trajNN knl', 'DoG knl', 'adj DoG knl'};
%}
algodirs={...
    'SE',...
    'nn',...
    'kernel',...
    'hwadjust'
    };
algofiles={...
    'originalbotrace_SEard_dim_50smooth',...
    'without_boombotrace_nn13_dim_50smooth',...
    'without_boombotrace_kernel_dim_50smooth_',...
    'without_boombotrace_kernel_dim_50smooth_hwadjust1_'
    };
clrs={'r', 'g', 'k', 'b'};
markers={'none', '^', '*', '+'};
lstyles={'-', '-', '-', '-'};
labels={'SE kernel', 'trajNN knl', 'DoG knl', 'adj DoG knl'};
%}

hs={};
nruns = 50;
n_trials = 50;
use_ci = true;

for di=1:length(algofiles)
    fname = char(strcat(dir_prefix, 'botrace_', algodirs{di}, ...
        '/atrias_', algofiles{di}, 'runs1_50.mat'));    
    fprintf('Loading data for %s\n', fname);
    assert(exist(fname, 'file') == 2)
    res = load(fname);
    Atranspose = res.botrace_50d.values;
    A = Atranspose';  % make rows be runs
    A = A(:,1:n_trials);  % clip to n_trials
    mn = min(A,[],2);
    A(mn==failed_cost,:)=NaN;
    nfailed_runs = sum(mn==failed_cost);
    if (nfailed_runs > 0)
        fprintf('failed runs: ')
        disp(find(mn==failed_cost)')
    end
    spec = {clrs{di}, markers{di}, lstyles{di}};
    [A,B,h] = plot_botrace_helper(NaN, spec, [cost_min, cost_max], use_ci, A);
    hs{di}=h;
    % Print various stats.
    %disp(A)
    tmp=B(:,50)
    disp(nanmean(tmp(tmp<cost_good)));
    disp(sum(tmp<cost_good)/sum(~isnan(tmp)));
end
fontsize=35;
l = legend([hs{1}, hs{2}, hs{3}, hs{4}], ...
    labels, 'FontSize',fontsize, 'Location','northeast','Interpreter','none');
xlabel('trials')
ylbl=ylabel('best cost so far');
set(ylbl, 'Units', 'Normalized', 'Position', [-0.06, 0.5, 0]);
xt = get(gca, 'XTick'); set(gca, 'FontSize', fontsize)
yt = get(gca, 'YTick'); set(gca, 'FontSize', fontsize)
set(gca,'fontsize',fontsize)
set(gcf,'color','w');

