% mkey = sprintf('%0.4f%0.4f%0.4f%0.4f%0.4f', z(i,1), z(i,2), z(i,3), z(i,4), z(i,5));
function [minsample,minvalue,botrace] = bayesoptBARE(...
    F,opt,rtp_or_model,rn,cost_op,num_rand_pts,scoresMap,walking_cost)
    if(nargin < 3)
        rtp_or_model = [];
        rn = [];
    end
    if (nargin < 6)
        num_rand_pts = 2;  % num random pts to request at the start of BO        
    end
    if(nargin < 7), scoresMap = containers.Map; end;  % empty map
    if(nargin < 8), walking_cost = -inf; end;  % empty map
    warning('on');
    % Check options for minimum level of validity
    check_opts(opt);

    % Either make a grid from a Sobol sequence or use opt.grid
    if isfield(opt,'grid')
        %hyper_grid = scale_point(opt.grid,opt.mins,opt.maxes);
        % Our grid is already scaled
        hyper_grid = opt.grid;
		opt.grid_size = size(hyper_grid,1);
    else
    	sobol = sobolset(opt.dims);
    	hyper_grid = sobol(1:opt.grid_size,:);
        % If the user wants to filter out some candidates
        if isfield(opt,'filter_func')
			hyper_grid = scale_point(opt.filter_func(...
                unscale_point(hyper_grid,opt.mins,opt.maxes)), ...
                opt.mins,opt.maxes);
        end
    end

    incomplete = logical(ones(size(hyper_grid,1),1));
    
    samples = [];
    values = [];
    hwadjustinfo = [];
    % tracking BO convergence
    minvalues = [];

    init = floor(rand(1,num_rand_pts)*opt.grid_size);
    fprintf('Running first point...\n');
    pt1 = unscale_point(hyper_grid(init(1),:),opt.mins,opt.maxes);
    disp_map(scoresMap, hyper_grid(init(1),:));

    stopTime = opt.stopTime;
    % Get values for the first two samples (no point using a GP yet)
    [value, ~, dog_score, ~, ~] = F(pt1, stopTime, rtp_or_model, cost_op);
    values = [values;value];
    if isfield(opt,'hwadjust')
        adjusted_dog_score = dog_adjust(dog_score, opt);
        hwadjustinfo = [hwadjustinfo;adjusted_dog_score]; 
    end;

    if num_rand_pts > 1
        assert(num_rand_pts == 2)
        fprintf('Running second point...\n');
        pt2 = unscale_point(hyper_grid(init(2),:),opt.mins,opt.maxes);
        disp_map(scoresMap, hyper_grid(init(2),:));
        [value, ~, dog_score, ~, ~] = F(pt2, stopTime, rtp_or_model, cost_op);
        values = [values;value];
        if isfield(opt,'hwadjust')
            adjusted_dog_score =  dog_adjust(dog_score, opt);
            hwadjustinfo = [hwadjustinfo; adjusted_dog_score];
        end
    end
    
    incomplete(init) = false;
    samples = [samples;hyper_grid(init,:)];
    
    % Remove first sample(s) from grid
    hyper_grid = hyper_grid(incomplete,:);
    if isfield(opt,'nn_out_vals')
        opt.nn_out_vals = opt.nn_out_vals(incomplete,:);
    end
    if isfield(opt,'nn_vars')
        opt.nn_vars = opt.nn_vars(incomplete,:);
    end
    incomplete = logical(ones(size(hyper_grid,1),1));

    % Main BO loop
	i_start = length(values) - num_rand_pts + 1;
    for i = i_start:opt.max_iters-num_rand_pts,hidx = -1;
	    % Score points on the grid using the acquisition function.
        [hyper_cand,hidx,aq_val] = get_next_cand(...
            samples,values,hyper_grid,opt,scoresMap,hwadjustinfo);
        
        disp(hyper_cand)
        if opt.use_ucb
            fprintf('Iteration %d, bnd = %f\n',i+num_rand_pts,-aq_val);
        else
            fprintf('Iteration %d, eic = %f\n',i+num_rand_pts,aq_val);
        end

        % Evaluate the candidate with the highest EI to get the actual function 
        % value, and add this function value and the candidate to our set.
        hyper_cand_scaled = scale_point(hyper_cand,opt.mins,opt.maxes);
        disp_map(scoresMap, hyper_cand_scaled(:));
        tic;
        [value, ~, dog_score, ~, ~] = F(hyper_cand, stopTime, rtp_or_model, cost_op);
        samples = [samples;scale_point(hyper_cand,opt.mins,opt.maxes)];
        values(end+1,:) = value;
        if isfield(opt,'hwadjust')
            adjusted_dog_score = dog_adjust(dog_score, opt);
            hwadjustinfo(end+1,:) = adjusted_dog_score;
        end

        % Remove this candidate from the grid (I use the incomplete vector like 
        % this because I will use this vector for other purposes in the future.)
        if hidx >= 0
            incomplete(hidx) = false;
        	hyper_grid = hyper_grid(incomplete,:);
            if isfield(opt,'nn_out_vals')
                opt.nn_out_vals = opt.nn_out_vals(incomplete,:);
            end
            if isfield(opt,'nn_vars')
                opt.nn_vars = opt.nn_vars(incomplete,:);
            end
        	incomplete = logical(ones(size(hyper_grid,1),1));
        end

        minvalues = [minvalues; min(values)];
	    fprintf(' value = %f, overall min = %f, hidx = %d\n',...
            value,min(values),hidx);
        botrace.samples = unscale_point(samples,opt.mins,opt.maxes);
        botrace.values = values;
        botrace.hwadjustinfo = hwadjustinfo;

        if opt.save_trace
            save(opt.trace_file,'botrace');
        end
        % For challenging settings it takes too long to simulat all the trials, 
        % so we only count the number of runs that found walking pts.
        if min(values) <= walking_cost, break; end;
    end
	
	% Get minvalue and minsample
    [mv,mi] = min(values);
    minvalue = mv;
    minsample = unscale_point(samples(mi,:),opt.mins,opt.maxes);
end


function [hyper_cand,hidx,aq_val] = get_next_cand(...
    samples,values,hyper_grid,opt,scoresMap,hwadjustinfo)
    % Get posterior means and variances for all points on the grid.
    if isfield(opt,'fixed_dims')
        assert(size(opt.fixed_dims,2) == opt.dims)
        % If opt.fixed_dims is [NaN, 3.0, NaN, 0.5] then will optimize only
        % over 1st and 3rd dimensions, then fill in 3.0 for 2nd, 0.5 for 4th.
        % Note: this could generate points off the grid.
        fixed_dim_ids = ~isnan(opt.fixed_dims);
        hyper_grid_backup_fixed = hyper_grid(:,fixed_dim_ids);
        scaled_fixed_dims = scale_point(opt.fixed_dims,opt.mins,opt.maxes);
        hyper_grid(:,fixed_dim_ids) = scaled_fixed_dims(fixed_dim_ids);
    end
    [mu,sigma2,hyp] = get_posterior(...
        samples,values,hyper_grid,opt,-1,scoresMap,hwadjustinfo);
    % Compute acquition function for all points in the grid and find maximum.
    best = min(values);
    ac = acquisition(best,mu,sigma2,opt);

	%[mei,meidx] = max(ac);
    if (isfield(opt,'random_topn') && opt.random_topn)  % avoid only getting one
        topn_n = min(1000, opt.grid_size);              % lucky point as maximum
        topn_ids = zeros(topn_n,1);                     % instead pick one of 1K
        tmp_ac = ac(:,1);                               % best points randomly
        for i=1:topn_n
            [~, topn_ids(i)] = max(tmp_ac);
            tmp_ac(topn_ids(i)) = -inf;
        end
        fprintf('random_topn:\n');
        disp(ac(topn_ids(1:10)));
        fprintf('mu:\n');
        disp(mu(topn_ids(1:10)));
        fprintf('sigma:\n');
        disp(sqrt(sigma2(topn_ids(1:10))));
        if isfield(opt, 'nn_out_vals')
            fprintf('nn_out_vals:\n');
            disp(opt.nn_out_vals(topn_ids(1:10),:))
        end
        if isfield(opt, 'nn_vars')
            fprintf('nn_vars:\n');
            disp(opt.nn_vars(topn_ids(1:10)))
        end
    else
        mei = max(ac);
        topn_ids = ac==mei;
    end
    aa = 1:length(ac);
    aa = aa(topn_ids);
    meidx = aa(randi([1,length(aa)]));
    fprintf('Number of best points %d best id %d \n', length(aa), meidx);
    hyper_cand = unscale_point(hyper_grid(meidx,:),opt.mins,opt.maxes);
    if isfield(opt,'fixed_dims')
        hyper_grid(:,fixed_dim_ids) = hyper_grid_backup_fixed(:,:);
    end
    hidx = meidx;
    aq_val = ac(meidx);
end

function [hyper_cand,hidx,aq_val] = get_next_cand_conditional(...
    samples,values,hyper_grid,cond,opt,scoresMap,hwadjustinfo)
    % Get next candidate conditioned on one dim fixed
    % cond is the last dimension - eg desired speed 

    idx = (hyper_grid(:,end) == cond);
    % Get posterior means and variances for all points on the grid that have
    % the last dimension same as cond.
    [mu,sigma2,hyp] = get_posterior(...
        samples,values,hyper_grid(idx,:),opt,-1,scoresMap,hwadjustinfo);

    % Compute EI for all points in the grid that have last dim as cond, 
    % and find the maximum.
    best = min(values);
    ac = acquisition(best,mu,sigma2,opt);

    [mei,meidx] = max(ac);
    hyper_cand = unscale_point(hyper_grid(meidx,:),opt.mins,opt.maxes);
    hidx = meidx;

    aq_val = mei;
end

function [mu,sigma2,hyp] = get_posterior(X,y,x_hats,opt,hyp,...
    scoresMap,hwadjustinfo)
    % Fit DoG hardware vs simluation adjustment GP.
    if length(hwadjustinfo) > 1
        fprintf('hwadjustinfo\n');
        disp(hwadjustinfo);
        assert(opt.hwadjust == 1)  % make sure we were asked to compute this
        fprintf('compute adjust_hyp\n');
        tic
        adjust_hyp = [];
        adjust_n_mh = num_hypers(@meanLocomotion,opt);  % needs simScoresMap
        adjust_n_ch = num_hypers(@covSEard,opt);
        adjust_hyp.mean = zeros(adjust_n_mh,1);
        adjust_hyp.cov = zeros(adjust_n_ch,1);  % cov will exp, so s^2 & length scales are 1.0
        adjust_hyp.lik = log(0.1);
        adjust_hyp = minimize(adjust_hyp,@gp,-100,opt.inffunc,...
            @meanLocomotion,@covSEard,@likGauss,X,hwadjustinfo);
        toc
        fprintf('optimized adjust_hyp: mean %0.4f lik %0.4f cov:\n', ...
            adjust_hyp.mean, adjust_hyp.lik);
        disp(adjust_hyp.cov);
        fprintf('getting adjust_mu\n');
        tic
        [adjust_mu,~] = gp(adjust_hyp,opt.inffunc,...
            @meanLocomotion,@covSEard,@likGauss,X,hwadjustinfo,opt.grid);
        toc
        fprintf('update scoresMaps with adjust_mu\n');
        tic
        % Update adjust estimates in scoresMap
        % Modify the mean by a 'deterministic' mean function composed of
        % DoG scores from simulation (see GPML book sec. 2.7 p.46).
        for i = 1:opt.grid_size
            mkey = sprintf('%0.4g', opt.grid(i,:));
            scoresMap(mkey) = adjust_mu(i);
        end
        toc
        fprintf('done updating scoresMaps\n');
    end
    % Get posterior mu and sigma2 for x_hats.
    if isfield(opt,'composite_gp')
        [mu,sigma2,hyp] = get_posterior_composite(...
            X,y,x_hats,opt,hyp,scoresMap);
    else
        [mu,sigma2,hyp] = get_posterior_gpml(X,y,x_hats,opt,hyp);
    end
end

function [mu,sigma2,hyp] = get_posterior_gpml(X,y,x_hats,opt,hyp)
    meanfunc = opt.meanfunc;
    covfunc = opt.covfunc;
    inffunc = opt.inffunc;
    if isnumeric(hyp)
        if isfield(opt,'num_mean_hypers')
            n_mh = opt.num_mean_hypers;
        else
            n_mh = num_hypers(meanfunc{1},opt);
        end
        if isfield(opt,'num_cov_hypers')
            n_ch = opt.num_cov_hypers;
        else
            n_ch = num_hypers(covfunc{1},opt);
        end
        hyp = [];
        hyp.mean = zeros(n_mh,1);
        hyp.cov = zeros(n_ch,1);  % cov will exp, so s^2 & length scales are 1.0
        hyp.lik = log(0.1);
        if size(y,1) > 1  % optimize hyp only after 2 points are obtained
            % Figure out whether we need to optimize hyper parameters.
            % If all the costs are the same - then it makes sense to keep the
            % default values, until more variable costs are obtained.
            % Value of std(y)=50 works for raibert 9D cost for highly
            % sample-efficient methods like DoG kernel.
            if (~isfield(opt,'hyp_min_std') || std(y) > 50)
                hyp = minimize(hyp,@gp,-100,inffunc,meanfunc,...
                    covfunc,@likGauss,X,y);
                fprintf('optimized hyp: mean %0.4f lik %0.4f cov:\n', ...
                    hyp.mean, hyp.lik);
                disp(hyp.cov);
            end
        end
    end
    [mu,sigma2] = gp(hyp,inffunc,meanfunc,covfunc,@likGauss,X,y,x_hats);
end

function ac = acquisition(best,mu,sigma2,opt)
    if (isfield(opt,'use_ucb') && opt.use_ucb)
        ac = compute_ucb(mu,sigma2);
    elseif (isfield(opt,'use_eri') && opt.use_eri)
        ac = compute_eri(best,mu,sigma2,opt.nn_vars);
    else
        ac = compute_ei(best,mu,sigma2);
    end
end

function ei = compute_ei(best,mu,sigma2)
    sigmas = sqrt(sigma2);
    u = (best - mu) ./ sigmas;
    ucdf = normcdf(u);
    updf = normpdf(u);
    ei = sigmas .* (u .* ucdf + updf);
end

function eri = compute_eri(best,mu,sigma2,risk_vars)
    % Expected Risk Improvement.
    % Adapted from "Variable Risk Control via Stochastic Optimization" by
    % Scott R. Kuindersma, Roderic A. Grupen, and Andrew G. Barto.
    k = 1.0;
    sigmas = sqrt(sigma2);
    u = (best - mu - k*risk_vars)./sigmas;
    ucdf = normcdf(u);
    updf = normpdf(u);
    eri = sigmas .* (u .* ucdf + updf);
    % Debug prints.
    fprintf('ERI\n');
    disp(u(1:10))
    orig_u = (best-mu)./sigmas;
    disp(orig_u(1:10));
    disp(risk_vars(1:10));
end

function bnd = compute_ucb(mu,sigma2)
    % From "Sequential Design Optimization" Cox and John 1992
    % select minimum over mu + k sigma
    % We are doing minimization, so technically will be using LCB
    % (lower confidence bound).
    % However, to make this consistent with EI, we will return the
    % negative of the bounds, hence the candidate with the lowest bound
    % would have the maximum acquisition value (ac).
    %
    k = 1.0; % Note: default from SheffieldML GPyOpt is 2; tested: 0.1-5.0
    sigmas = sqrt(sigma2);
    bnd = -1.0*(mu - k*sigmas);
end

function bnd = compute_ucb_risk_averse(mu,sigma2)
    % A simple Confidence Bound Criterion for risk-averse experiments.
    % See "Variable Risk Control via Stochastic Optimization" from 2013 by
    % Scott R. Kuindersma, Roderic A. Grupen, and Andrew G. Barto.
    k = 1.0;
    sigmas = sqrt(sigma2);
    bnd = -1.0*(mu + k*sigmas);
end

function upt = unscale_point(x,mins,maxes)
    if size(x,1) == 1
        upt = x .* (maxes - mins) + mins;
    else
        upt = bsxfun(@plus,bsxfun(@times,x,(maxes-mins)),mins);
    end
end

function pt = scale_point(x,mins,maxes)
	pt = bsxfun(@rdivide,bsxfun(@minus,x,mins),maxes-mins);
end
    
function check_opts(opt)
    if ~isfield(opt,'dims')
        error('bayesopt:opterror',['The dims option specifying the dimensionality of ' ...
            'the optimization problem is required']);
    end
    
    if ~isfield(opt,'mins') || length(opt.mins) < opt.dims
        error('bayesopt:opterror','Must specify minimum values for each hyperparameter');
    end
       
    if ~isfield(opt,'maxes') || length(opt.maxes) < opt.dims
        error('bayesopt:opterror','Must specify maximum values for each hyperparameter');
    end 
end


    
function nh = num_hypers(func,opt)
    str = func();
    nm = str2num(str);
    if ~isempty(nm)
        nh = nm;
    else
        if all(str == '(D+1)')
            nh = opt.dims + 1;
        elseif all(str == '(D+2)')
            nh = opt.dims + 2;
        else
            error('bayesopt:unkhyp','Unknown number of hyperparameters asked for by one of the functions');
        end
    end
end

function [res] = disp_map(scoresMap, x)
    if isa(scoresMap, 'py.nnmod.NNMod')
        nn_outs = scoresMap.get_nn_outs(x);
        nn_outs_m = double(py.array.array('d',py.numpy.nditer(nn_outs)));
        fprintf('TF NN info: nnvar %0.4f nnouts ', nn_outs_m(:,end))
        res = nn_outs_m(:,1:end-1);
        disp(res);
    elseif ~scoresMap.isempty()
        mkey = sprintf('%0.4g', x);
        assert(scoresMap.isKey(mkey));
        fprintf('Map info');
        res = scoresMap(mkey);
        disp(res);
    else
        res = NaN;
    end
end

% Adjust dog_score computed from a hardware-like run to be comparable to
% simulation-based score obatined from 5s runs.
function [adjusted_dog_score] = dog_adjust(hw_dog_score, opt)
    if ~isfield('sim_stopTime',opt)
        adjust_ratio = 5/opt.stopTime;  % usually 5s simulations, so 5/30=1/6
    else
        adjust_ratio = opt.sim_stopTime/opt.stopTime;
    end
    adjusted_dog_score = hw_dog_score*adjust_ratio;
end
