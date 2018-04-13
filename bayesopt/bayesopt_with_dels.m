function [minsample,minvalue,botrace] = bayesopt_with_dels(F,opt,dels, vis)

% TODO(rika): comment out when debugging
    % Seeing warnings from chol. Maybe do something to get rid of these.
    % They cause unstable estimates of the acquisition function.
    warning('on');
    % Check options for minimum level of validity
    check_opts(opt);

    % Draw initial candidate grid from a Sobol sequence
	if isfield(opt,'grid')
		%hyper_grid = scale_point(opt.grid,opt.mins,opt.maxes);
        % Our grid is already scaled
        hyper_grid = opt.grid;
		opt.grid_size = size(hyper_grid,1);
	else 
    	sobol = sobolset(opt.dims);
    	hyper_grid = sobol(1:opt.grid_size,:);
		if isfield(opt,'filter_func'), % If the user wants to filter out
			hyper_grid = scale_point(opt.filter_func(...  % some candidates
                unscale_point(hyper_grid,opt.mins,opt.maxes)),opt.mins,opt.maxes);
		end
    end
    if(nargin < 3)
        dels = zeros(opt.grid_size,1);
    end
    if(nargin < 4)
        vis = 0;
    end
    
    incomplete = logical(ones(size(hyper_grid,1),1));
    hyper_grid_old = hyper_grid;
    
    samples = [];
    values = [];
    times = [];
    % tracking BO convergence
    minvalues = [];

    init = floor(rand(1,2)*opt.grid_size);
    fprintf('Running first point...\n');
    pt1 = unscale_point(hyper_grid(init(1),:),opt.mins,opt.maxes);
    % Get values for the first two samples (no point using a GP yet)
    [val1] = F(pt1);

    fprintf('Running second point...\n');
    pt2 = unscale_point(hyper_grid(init(2),:),opt.mins,opt.maxes);
    [val2] = F(pt2);
    
    incomplete(init) = false;
    samples = [samples;hyper_grid(init,:)];
    values = [values;val1;val2];
    
    if(vis)
        % tracking how BO is working
        subplot(3,1,1);
        plot(pt1,val1,'g*');
        hold on
        plot(pt2,val2,'*g');
        drawnow;
    end
    
    % Remove first two samples from grid
    hyper_grid = hyper_grid(incomplete,:);
    %incomplete = logical(ones(size(hyper_grid,1),1));
    
    % Main BO loop
	i_start = length(values) - 2 + 1;
    ids = 1:opt.grid_size;
    for i = i_start:opt.max_iters-2,hidx = -1;
	    % Score points on the grid using the acquisition function.
        idz = ids(incomplete);
        idx = ids(~incomplete);
        [hyper_cand,hidx,aq_val] = get_next_cand(samples,values,hyper_grid,opt, dels,idx,idz, vis);

        if opt.use_ucb
            fprintf('Iteration %d, bnd = %f',i+2,-aq_val);
        else
            fprintf('Iteration %d, eic = %f',i+2,aq_val);
        end

        % Evaluate the candidate with the highest EI to get the actual function 
        % value, and add this function value and the candidate to our set.
        tic;
        value = F(hyper_cand);
        times(end+1) = toc;
        samples = [samples;scale_point(hyper_cand,opt.mins,opt.maxes)];
        values(end+1) = value;
        
        		
        if(vis)
            % tracking how BO is working
            subplot(3,1,1);
            plot(hyper_cand,values(end),'g*');
            drawnow;
        end

        % Remove this candidate from the grid (I use the incomplete vector like 
        % this because I will use this vector for other purposes in the future.)
        if hidx >= 0
            tmp = incomplete(incomplete);
        	tmp(hidx) = false;
            incomplete(incomplete) = tmp;
        	hyper_grid = hyper_grid_old(incomplete,:);
        	% incomplete = logical(ones(size(hyper_grid,1),1));
        end
		
        minvalues = [minvalues; min(values)];
	    fprintf(', value = %f, overall min = %f, hidx = %d\n',value,min(values),hidx);
        botrace.samples = unscale_point(samples,opt.mins,opt.maxes);
        botrace.values = values;
        botrace.times = times;
        
        if(vis)
            % track BO convergence
            subplot(3,1,3);
            plot(minvalues,'r*-');
            title('fn min value');
            drawnow;
        end

        if opt.save_trace
            save(opt.trace_file,'botrace');
        end
    end
	
	% Get minvalue and minsample
    [mv,mi] = min(values);
    minvalue = mv;
    minsample = unscale_point(samples(mi,:),opt.mins,opt.maxes);
end


function [hyper_cand,hidx,aq_val] = get_next_cand(samples,values,hyper_grid,opt,dels,idx,idz,vis)
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
    [mu,sigma2,ei_hyp] = get_posterior(samples,values,hyper_grid,opt,-1, dels,idx,idz);
    % Compute acquition function for all points in the grid and find maximum.
    best = min(values);
    ac = acquisition(best,mu,sigma2,opt.use_ucb);
    
    if(vis)
        subplot(3,1,2); plot(unscale_point(hyper_grid,opt.mins,opt.maxes),ac,'.');
        drawnow;
        f = [mu+2*sqrt(sigma2); flip(mu-2*sqrt(sigma2),1)];
        xs = linspace(opt.mins, opt.maxes, opt.grid_size)';
        subplot(3,1,1); plot(xs(idz), mu)
        fill([xs(idz); flip(xs(idz),1)], f, [7 7 7]/8)
    end
    
    [mei,meidx] = max(ac);
    hyper_cand = unscale_point(hyper_grid(meidx,:),opt.mins,opt.maxes);
    if isfield(opt,'fixed_dims')
        hyper_grid(:,fixed_dim_ids) = hyper_grid_backup_fixed(:,:);
    end
    hidx = meidx;
    aq_val = mei;
end

function [hyper_cand,hidx,aq_val] = get_next_cand_conditional(...
    samples,values,hyper_grid, cond, opt)
    % Get next candidate conditioned on one dim fixed
    % cond is the last dimension - eg desired speed 

    idx = (hyper_grid(:,end) == cond);
    % Get posterior means and variances for all points on the grid that have
    % the last dimension same as cond.
    [mu,sigma2,ei_hyp] = get_posterior(samples,values,hyper_grid(idx,:),opt,-1);

    % Compute EI for all points in the grid that have last dim as cond, 
    % and find the maximum.
    best = min(values);
    ac = acquisition(best,mu,sigma2,opt.use_ucb);

    [mei,meidx] = max(ac);
    hyper_cand = unscale_point(hyper_grid(meidx,:),opt.mins,opt.maxes);
    hidx = meidx;

    aq_val = mei;
end

function [mu,sigma2,hyp] = get_posterior(X,y,x_hats,opt,hyp,dels,idx,idz)
    %Akshara : adding noise for simulation data
    global delx delz
    delx = dels(idx)';
    delz = dels(idz);
    
    meanfunc = opt.meanfunc;
    covfunc = opt.covfunc;
    if isnumeric(hyp)
        if isfield(opt,'num_mean_hypers'),
            n_mh = opt.num_mean_hypers;
        else
            n_mh = num_hypers(meanfunc{1},opt);
        end
        if isfield(opt,'num_cov_hypers'),
            n_ch = opt.num_cov_hypers;
        else
            n_ch = num_hypers(covfunc{1},opt);
        end
        hyp = [];
        hyp.mean = zeros(n_mh,1);
        hyp.cov = zeros(n_ch,1);
        hyp.lik = log(0.1);
		hyp = minimize(hyp,@gp,-100,@infExact,meanfunc,covfunc,@likGauss,X,y);
    end
    [mu,sigma2] = gp(hyp,@infExact,meanfunc,covfunc,@likGauss,X,y,x_hats);
%     fprintf('INSIDE: get_posterior() hyp.mean=%0.2f hyp.lik=%0.2f hyp.cov=\n', ...
%         hyp.mean, hyp.lik);
%     fprintf(' %0.2f', hyp.cov);
%     fprintf('\n');
end

function ac = acquisition(best,mu,sigma2,use_ucb)
    if use_ucb
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
    k = 1.0; % TODO(rika) default from SheffieldML GPyOpt is 2; tested: 0.1-5.0
    sigmas = sqrt(sigma2);
    bnd = -mu + k*sigmas;
end

function upt = unscale_point(x,mins,maxes)
    if size(x,1) == 1,
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


