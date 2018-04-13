
%% Analyze explicit mismatch simulation experiments.

% Plots for BO results with error region.
addpath(genpath('../../tools'))  % for plot_botrace_helper

%% Load BO runs data.
dir_prefix = '../../../boforloco_runs_writeups/runs_mismatch/';
%cost_max = 1.3; cost_good = 0.6; failed_cost = NaN;
cost_max = 100; cost_good = 20; failed_cost = 120;

extra_clrs = struct(...
    'lo', [1, 0.85, 0.7], ...  % extra light orange
    'do', [1, 0.75, 0.2], ...  % light orange
    'edo', [1, 0.5, 0],...  % extra dark orange
    'elm', [1.0, 0.7, 1.0], ...  % extra light magenta
    'lm', [0.9, 0, 0.9], ...  % light magenta
    'dm', [0.4, 0, 0.4],...  % extra dark magenta
    'elb', [0.7, 0.7, 1.0], ...  % extra light blue
    'lb', [0, 0, 0.9], ...  % light blue
    'db', [0, 0, 0.4],...  % extra dark blue
    'ly', [0.9, 0.9, 0], ...  % light yellow
    'dy', [0.5, 0.5, 0], ...  % dark yellow
    'edy', [0.3, 0.3, 0]);  % extra dark yellow
%{
algodirs={...
    'prior',...
    'prior',...
    'prior',...
    'SE'
    };
algofiles={...
    'without_boombotrace_prior_dim_50',...
    'gear_dynamicsbotrace_prior_dim_50',...
    'originalbotrace_prior_dim_50',...
    'originalbotrace_SEard_50'
    };
clrs={extra_clrs.edo, extra_clrs.do, extra_clrs.lo, 'r'};
markers={'>', 'p', 'none' 'none'};
lstyles={'-', '-', '--', '-'};
labels={['Cost prior', newline, 'built using:' newline, ...
    'no boom'], ...
    'simple gears', ...
    'original sim.'};
labels2={'SE kernel'};
%}
%{
algodirs={...
    'cully',...
    'cully',...
    'cully',...
    'SE'
    };
algofiles={...
    'without_boombotrace_cully_prior_dim_50_',...
    'gear_dynamicsbotrace_cully_prior_dim_50_',...
    'originalbotrace_cully_prior_dim_50_',...
    'originalbotrace_SEard_50'
    };
clrs={extra_clrs.edy, extra_clrs.dy, extra_clrs.ly, 'r'};
markers={'>', 'p', 'none' 'none'};
lstyles={'-', '-', '--', '-'};
labels={['M-BOA (w/ prior)', newline, 'built using:' newline, ...
    'no boom'], ...
    'simple gears', ...
    'original sim.'};
labels2={'SE kernel'};
%}
%{
algodirs={...
    'prior_and_kernel',...
    'prior_and_kernel',...
    'prior_and_kernel',...
    'SE'
    };
algofiles={...
    'without_boombotrace_prior_and_kernel_50hwadjust1_',...
    'gear_dynamicsbotrace_prior_and_kernel_dim_50hwadjust1_',...
    'originalbotrace_prior_and_kernel_dim_50hwadjust1_',...
    'originalbotrace_SEard_50'
    };
clrs={'k>-', 'bp-', 'g--', 'r-'};
labels={['adj DoG knl w/ 5s cost prior using :' newline, ...
    'no boom & simple gears'], ...
    'boom & simple gears', ...
    'original simulator'};
labels2={'SE kernel'};
%}
%{
algodirs={...
    'cully',...
    'cully',...
    'cully',...
    'SE'
    };
algofiles={...
    'without_boombotrace_cully_dim_50_',...
    'gear_dynamicsbotrace_cully_dim_50_',...
    'originalbotrace_cully_dim_50_',...
    'originalbotrace_SEard_50'
    };
clrs={extra_clrs.dm, extra_clrs.lm, extra_clrs.elm, 'r'};
markers={'>', 'p', 'none' 'none'};
lstyles={'-', '-', '--', '-'};
labels={['M-BOA-knl', newline, 'built using:' newline, ...
    'no boom'], ...
    'simple gears', ...
    'original sim.'};
labels2={'SE kernel'};
%}

algodirs={...
    'hwadjust',...
    'hwadjust',...
    'hwadjust',...
    'SE'
    };
algofiles={...
    'without_boombotrace_kernel_dim_50hwadjust1_',...
    'gear_dynamicsbotrace_kernel_dim_50hwadjust1_',...
    'originalbotrace_kernel_dim_50hwadjust1_',...
    'originalbotrace_SEard_50'
    };
clrs={extra_clrs.db, extra_clrs.lb, extra_clrs.elb, 'r'};
markers={'>', 'p', 'none' 'none'};
lstyles={'-', '-', '--', '-'};
labels={['adj DoG knl', newline, 'built using:' newline, ...
    'no boom'], ...
    'simple gears', ...
    'original sim.'};
labels2={'SE kernel'};


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
    if strcmp(algodirs{di}, 'hyp_default'), res.botrace_50d = res.botrace_5d; end
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
    [A,B,h] = plot_botrace_helper(NaN, spec, cost_max, use_ci, A);
    hs{di}=h;
    % Print various stats.
    %disp(A)
    tmp=B(:,40);
    disp(nanmean(tmp(tmp<50)));
    disp(sum(tmp<50)/sum(~isnan(tmp)));
end
fontsize=35;
l = legend([hs{1}, hs{2}, hs{3}], ...
    labels, 'FontSize',fontsize, 'Location','northeast','Interpreter','none');
%l = legend([hs{1}, hs{2}, hs{3}, hs{4}], ...
%    labels, 'FontSize',fontsize, 'Location','northeast','Interpreter','none');
xlabel('trials')
ylbl=ylabel('best cost so far');
set(ylbl, 'Units', 'Normalized', 'Position', [-0.06, 0.5, 0]);
xt = get(gca, 'XTick'); set(gca, 'FontSize', fontsize)
yt = get(gca, 'YTick'); set(gca, 'FontSize', fontsize)
set(gca,'fontsize',fontsize)
set(gcf,'color','w');
% Second legend.
% https://se.mathworks.com/matlabcentral/answers/365006-how-to-create-2-legends-in-one-figure
a=axes('position',get(gca,'position'),'visible','off');
legend(a,[hs{4}],labels2,'FontSize',fontsize,'Location','northeast');

