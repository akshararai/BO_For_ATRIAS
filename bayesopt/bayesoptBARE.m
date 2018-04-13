% mkey = sprintf('%0.4f%0.4f%0.4f%0.4f%0.4f', z(i,1), z(i,2), z(i,3), z(i,4), z(i,5));
function [minsample,minvalue,botrace] = bayesoptBARE(...
    F,opt,rtp_or_model,cost_op,num_rand_pts,scoresMapF,siminit_name)
    %warning('on');
    check_opts(opt);  % check options for minimum level of validity
    if isfield('hwadjust', opt)
        assert(isfield('hwadjust_orig_sim_dog_scores', opt));
    end
    % Initialize locomotion-specific arguments.
    if ~isfield('stopTime',opt)
        stopTime = 30;
    else
        stopTime = opt.stopTime;
    end
    if ~isfield('kernel_dims',opt)
        opt.kernel_dims = opt.dims;
    end
    if (nargin < 3)
        rtp_or_model = [];
        rn = [];
    end
    if (nargin < 5), cost_op = str2func('findCost_metabolic'); end
    if (nargin < 6)
        num_rand_pts = 2;  % num random pts to request at the start of BO        
    end
    if(nargin < 7), scoresMapF = @(x) NaN; end  % dummy map
    if(nargin < 8), siminit_name = []; end

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

    init = randi(opt.grid_size, num_rand_pts, 1)%floor(rand(1,num_rand_pts)*opt.grid_size);
    fprintf('Running first point...\n');
    pt1 = unscale_point(hyper_grid(init(1),:),opt.mins,opt.maxes);
    sim_score = disp_map(scoresMapF, hyper_grid(init(1),:));

    % Get values for the first two samples (no point using a GP yet)
    [value, ~, dog_score, ~, ~] = F(pt1, stopTime, rtp_or_model, cost_op, siminit_name);
    values = [values;value];
    if isfield(opt,'hwadjust')
        adjusted_dog_score =  dog_adjust(dog_score, opt, stopTime);
        hwadjustinfo = [hwadjustinfo;(sim_score(1)-adjusted_dog_score)/opt.hwadjust_scale]; 
    end

    if num_rand_pts > 1
        assert(num_rand_pts == 2)
        fprintf('Running second point...\n');
        pt2 = unscale_point(hyper_grid(init(2),:),opt.mins,opt.maxes);
        sim_score = disp_map(scoresMapF, hyper_grid(init(2),:));
        [value, ~, dog_score, ~, ~] = F(pt2, stopTime, rtp_or_model, cost_op, siminit_name);
        values = [values;value];
        if isfield(opt,'hwadjust')
            adjusted_dog_score =  dog_adjust(dog_score, opt, stopTime);
            hwadjustinfo = [hwadjustinfo; (sim_score(1)-adjusted_dog_score)/opt.hwadjust_scale];
        end
    end
    
    incomplete(init) = false;
    samples = [samples;hyper_grid(init,:)];
    
    % Remove first sample(s) from grid
    hyper_grid = hyper_grid(incomplete,:);
    incomplete = logical(ones(size(hyper_grid,1),1));

    % Main BO loop
	i_start = length(values) - num_rand_pts + 1;
    for i = i_start:opt.max_iters-num_rand_pts,hidx = -1;
	    % Score points on the grid using the acquisition function.
        [hyper_cand,hidx,aq_val] = get_next_cand(...
            samples,values,hyper_grid,opt,scoresMapF,hwadjustinfo);
        
%         disp(hyper_cand)
        if opt.use_ucb
            fprintf('Iteration %d, bnd = %f\n',i+num_rand_pts,-aq_val);
        else
            fprintf('Iteration %d, eic = %f\n',i+num_rand_pts,aq_val);
        end

        % Evaluate the ca.ndidate with the highest EI to get the actual function 
        % value, and add this function value and the candidate to our set.
        hyper_cand_scaled = scale_point(hyper_cand,opt.mins,opt.maxes);
        sim_score = disp_map(scoresMapF, hyper_cand_scaled(:)');
        tic;
        [value, ~, dog_score, ~, ~] = F(hyper_cand, stopTime, rtp_or_model, cost_op, siminit_name);
        samples = [samples;scale_point(hyper_cand,opt.mins,opt.maxes)];
        values(end+1,:) = value;
        if isfield(opt,'hwadjust')
            adjusted_dog_score = dog_adjust(dog_score, opt, stopTime);
            hwadjustinfo(end+1,:) = (sim_score(1)-adjusted_dog_score)/opt.hwadjust_scale;
        end

        % Remove this candidate from the grid (I use the incomplete vector like 
        % this because I will use this vector for other purposes in the future.)
        if hidx >= 0
            incomplete(hidx) = false;
        	hyper_grid = hyper_grid(incomplete,:);
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
    end
	
	% Get minvalue and minsample
    [mv,mi] = min(values);
    minvalue = mv;
    minsample = unscale_point(samples(mi,:),opt.mins,opt.maxes);
end


function [hyper_cand,hidx,aq_val] = get_next_cand(...
    samples,values,hyper_grid,opt,scoresMapF,hwadjustinfo)
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
        samples,values,hyper_grid,opt,-1,scoresMapF,hwadjustinfo);
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
    samples,values,hyper_grid,cond,opt,scoresMapF,hwadjustinfo)
    % Get next candidate conditioned on one dim fixed
    % cond is the last dimension - eg desired speed 

    idx = (hyper_grid(:,end) == cond);
    % Get posterior means and variances for all points on the grid that have
    % the last dimension same as cond.
    [mu,sigma2,hyp] = get_posterior(...
        samples,values,hyper_grid(idx,:),opt,-1,scoresMapF,hwadjustinfo);

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
    scoresMapF,hwadjustinfo)
    % Fit DoG hardware vs simluation adjustment GP.
    % For low-dimensional kernels like DoG1 we only want to enable computing
    % mismatch once some signal in the cost is observed. Otherwise we do not
    % know which regions of the space are better than others - cost of ~100
    % everywhere gives no signal. If the cost is ~100 everywhere, then the
    % best we can do for 1D kernel, for example, is to randomly sample in
    % the space that was already remapped to stretch the part of the space with
    % points that have a hope of walking and shrink the part of the space 
    % with non-walking points.
    has_signal = true;
    if isfield(opt, 'hwadjust_signal_threshold')
        has_signal = (max(y) - min(y)) > opt.hwadjust_signal_threshold;
    end
    if (length(hwadjustinfo) > 1)
        fprintf('hwadjustinfo\n');
        disp(hwadjustinfo*opt.hwadjust_scale);
        assert(opt.hwadjust >= 1)  % make sure we were asked to compute this
        fprintf('compute adjust_hyp\n');
        tic
        adjust_hyp = [];
        adjust_n_mh = num_hypers({@meanConst},opt.dims);
        adjust_n_ch = num_hypers({@covSEard},opt.dims);
        adjust_hyp.mean = zeros(adjust_n_mh,1);
        adjust_hyp.cov = zeros(adjust_n_ch,1);  % cov will exp, so s^2 & length scales are 1.0
        adjust_hyp.lik = log(0.1);
        if has_signal
            fprintf('optimizing adjust_hyp...\n')
            adjust_hyp = minimize(adjust_hyp,@gp,-100,opt.inffunc,...
            @meanConst,@covSEard,@likGauss,X,hwadjustinfo);
        end
        toc
        fprintf('adjust_hyp: mean %0.4f lik %0.4f cov:\n', ...
            adjust_hyp.mean, adjust_hyp.lik);
%         disp(adjust_hyp.cov')
        fprintf('getting adjust_mu\n');
        tic
        [adjust_mu,~] = gp(adjust_hyp,opt.inffunc,...
            @meanConst,@covSEard,@likGauss,X,hwadjustinfo,opt.grid);
        toc
        fprintf('update scoresMapF with adjust_mu\n');
        tic
        % Update adjust estimates in scoresMap
        res = [opt.hwadjust_orig_sim_dog_scores, adjust_mu];
        size(res)
        scoresMapF(opt.grid, res);
        toc
        fprintf('done updating scoresMapF\n');
    end
    [mu,sigma2,hyp] = get_posterior_gpml(X,y,x_hats,opt,hyp);
end

function [mu,sigma2,hyp] = get_posterior_gpml(X,y,x_hats,opt,hyp)
    meanfunc = opt.meanfunc;
    covfunc = opt.covfunc;
    inffunc = opt.inffunc;
    if isnumeric(hyp)
        if isfield(opt,'num_mean_hypers')
            n_mh = opt.num_mean_hypers;
        else
            n_mh = num_hypers(meanfunc,opt.kernel_dims);
        end
        if isfield(opt,'num_cov_hypers')
            n_ch = opt.num_cov_hypers;
        else
            n_ch = num_hypers(covfunc,opt.kernel_dims);
        end
        hyp = [];
        hyp.mean = zeros(n_mh,1);
        hyp.cov = zeros(n_ch,1);  % cov will exp, so s^2 & length scales are 1.0
        hyp.lik = log(0.1);
        if size(y,1) > 1  % optimize hyp only after 2 points are obtained
            has_signal = true;
            if isfield(opt, 'hwadjust_signal_threshold')
                has_signal = (max(y) - min(y)) > opt.hwadjust_signal_threshold;
            end
            if has_signal
                fprintf('optimizing hyp...\n')
                hyp = minimize(hyp,@gp,-100,inffunc,meanfunc,covfunc,@likGauss,X,y);
            end
            fprintf('hyp: mean %0.4f lik %0.4f cov:\n', hyp.mean, hyp.lik);
            disp(hyp.cov');
        end
    end
    [mu,sigma2] = gp(hyp,inffunc,meanfunc,covfunc,@likGauss,X,y,x_hats);
end

function ac = acquisition(best,mu,sigma2,opt)
    if (isfield(opt,'use_ucb') && opt.use_ucb)
        ac = compute_ucb(mu,sigma2);
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
    
function nh = num_hypers(funcstr,dims)
    str = feval(funcstr{:});
    if strcmp(str, '(D+1)')
        nh = dims + 1;
    elseif strcmp(str, '(D+2)')
        nh = dims + 2;
    else
        nh = eval(str);
    end
end

function [res] = disp_map(scoresMapF, x)
if isempty(scoresMapF) res=0; return; end 
    res = scoresMapF(x);
    fprintf('scoresMapF: ')
%     disp(x)
%     fprintf('->')
    disp(res)
end

% Adjust dog_score computed from a hardware-like run to be comparable to
% simulation-based score obatined from 5s runs.
function [adjusted_dog_score] = dog_adjust(hw_dog_score, opt, stopTime)
    if ~isfield('sim_stopTime',opt)
        sim_stopTime = 5;  % usually 5s simulations, so 5/30=1/6
    else
        sim_stopTime = opt.sim_stopTime;
    end
    adjust_ratio = sim_stopTime/stopTime;
    adjusted_dog_score = hw_dog_score*adjust_ratio;
end
