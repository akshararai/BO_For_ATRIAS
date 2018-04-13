% stopTime is timeout in seconds
function [cost1, cost2, total_score, score_4d, simout_vec] = evaluateGaitScore(...
    simout, stopTime, cost_op)

    if(nargin < 3)
        cost_op = [];
    end
    
    if(nargin <2)
        stopTime = 100;
    end
    
    cost1 = nan;
    cost2 = nan;
    total_score = 0;
    score_4d = zeros(4,1);
    simout_vec = [];

%    fprintf('evaluateGaitScore(): stopTime=%d sec\n', stopTime);

%% Record trajectory features from simout into a small (8d) simout_vec.
% simout_to_vec() wraps a try/catch, so safe to call on any input.
simout_vec = simout_to_vec(simout);

%% Scores

% VNMC_on = get(simout, 'VNMC_on');
% vnmc_tmp = VNMC_on.signals.values;
% vnmc(:) = vnmc_tmp(1,1,:);
% idxs = find(vnmc);
% if(length(idxs) < 1)
%     return;
% end
start_idx = 1;%idxs(1);
% diff = idxs(end) - idxs(1);

swingFootz = get(simout, 'PrimaryToSecondaryZ');
swingFootz = swingFootz.signals.values;
swingFootz = swingFootz(start_idx:end,:);

CoMz = get(simout, 'CoMHeight');
CoMz = CoMz.signals.values;
CoMz = CoMz(start_idx:end,:);

theta = get(simout, 'TorsoPitch');
theta = theta.signals.values;
theta = theta(start_idx:end,:);

CoMdx = get(simout, 'CoMVelocity');
CoMdx = CoMdx.signals.values;
CoMdx = CoMdx(start_idx:end,:);

% stanceFlag = get(simout, 'EnteringStance');
% stanceFlag = stanceFlag.signals.values;
LStance = get(simout, 'LStance');
LStance = LStance.signals.values;
LStance = LStance(start_idx:end,:);

RStance = get(simout, 'RStance');
RStance = RStance.signals.values;
RStance = RStance(start_idx:end,:);

rsteps = find(RStance);
lsteps = find(LStance);

if (length(LStance) < 30 || length(RStance) < 30)
    return;
end
if(length(rsteps) < 1 || length(lsteps) < 1)
    return;
end

step_idx(1) = min(rsteps(1), lsteps(1));
steps = find(LStance & RStance);
if(length(steps) < 1)
    return;
end

step_idx(2) = steps(1);
count = 3;
for i = 1:length(steps)-1
    if(steps(i+1) - steps(i) > 1)
        step_idx(count) = steps(i+1);
        count = count+1;
    end
end

nSteps = count -1;
v_avg = zeros(nSteps-1,1);
for step = 1:nSteps-1
    v_avg(step) = mean(CoMdx(step_idx(step):step_idx(step+1)));
end


if(nSteps >= 1)
%     Lflight = Rstep;
%     Rflight = Lstep;
    
    % knee flexion == ankle trajectory
    knee_score = 0;
    com_score = 0;
    hip_score = 0;
    for i=1:nSteps-1
        step_len = step_idx(i+1) - step_idx(i);
        if(step_len < 50)
            continue;
        end
        footz = swingFootz(step_idx(i):step_idx(i+1));
        comz = CoMz(step_idx(i):step_idx(i+1));
        pitch = theta(step_idx(i):step_idx(i+1));
        if(max(footz) > 0.1)
            knee_score = knee_score + 1;
        end
        if abs(mean(comz(1:20)) - mean(comz(end-20:end))) < 0.025
            com_score = com_score + 1;
        end
        if abs(mean(pitch(1:20)) - mean(pitch(end-20:end))) < 0.01
            hip_score = hip_score + 1;
        end
    end
%     knee_score = knee_score/2;
%     com_score = 1 - com_score;
%     hip_score = 0.5 - hip_score;
    
    % Mean velocity
    time = get(simout,'time');
    time = time.signals.values;

%     comx_score = v_avg(end);
    comx_score = sum(v_avg);
    total_score = knee_score + com_score + hip_score + comx_score; 
    total_score = total_score*time(end)/stopTime;
    score_4d = [knee_score, com_score, hip_score, comx_score];
    score_4d = score_4d.*time(end)/stopTime;
end


%% Costs
[cost1] = findCost(simout, stopTime, str2func('get'));
cost2 = findCost_metabolic(simout, stopTime);

if(~isempty(cost_op))
    cost1 = cost_op(simout, stopTime);
end

end