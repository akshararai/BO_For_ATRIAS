%% Test covWarp from gpml_2017 for DoG and NN kernel

addpath(genpath('../../bayesopt'));
rmpath(genpath('../../bayesopt/gpml'));
addpath(genpath('../../bayesopt/private/gpml_2017'));

%%
miter = -100;  %  minimize: if negative, |miter| is max allowed fxn evals
approx = 0;
warp = 1;

[x, y, xs] = gen_data1();

%meanfunc = [];                    % empty: don't use a mean function
meanfunc = {@meanSum, {@meanLinear, @meanConst}};  % affine mean fxn
if warp
  D = 1;  % dimensionality of the transformed input
  ddummy = @(x) 1;  % dummy derivative function
  grid = [x; xs];
  warped = nan(size(grid));
  for i=1:size(grid,1), warped(i,:) = locoTransformOne(grid(i,:)); end
  locomap([], [], grid, warped)
  covfunc = {@covWarp,{@covSEard},@locomap,ddummy,D};
else
  covfunc = @covSEard;
end
likfunc = @likGauss;              % Gaussian likelihood
hyp = struct('mean', [0 0], 'cov', [0 0], 'lik', log(0.1));

if approx
  u = linspace(-0.9,0.9,7)';                             % 1D grid of inducing pts
  %iu = randperm(size(x,1)); iu = iu(1:10); u = x(iu,:);  % subsample x
  hyp.xu = u;                                             % optimize u
  covfunc = {@apxSparse, {covfunc}, u};
  infr = @(varargin) infGaussLik(varargin{:}, struct('s', 1.0));
else
  u = [];
  infr = @infGaussLik;
end
hyp2 = minimize(hyp, @gp, miter, infr, meanfunc, covfunc, likfunc, x, y)
nlml = gp(hyp2, infr, meanfunc, covfunc, likfunc, x, y)
[mu, s2] = gp(hyp2, infr, meanfunc, covfunc, likfunc, x, y, xs);

plot_gp(mu, s2, x, y, xs, u)

%%

function [x, y, xs] = gen_data0()
  % Generate data from an example GP.
  % Second example from GPML tutorial:
  % http://www.gaussianprocess.org/gpml/code/matlab/doc/
  meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];
  covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
  hyp.lik = log(0.1);  % std of noise is 0.1

  n = 20;
  x = gpml_randn(0.3, n, 1);
  K = feval(covfunc{:}, hyp.cov, x);
  mu = feval(meanfunc{:}, hyp.mean, x);
  y = chol(K)'*gpml_randn(0.15, n, 1) + mu + exp(hyp.lik)*gpml_randn(0.2, n, 1);

  xs = linspace(-1.9, 1.9, 101)';  
  %plot(x, y, '+')
end


function [x, y, xs] = gen_data1()
  % First example from GPML tutorial:
  % http://www.gaussianprocess.org/gpml/code/matlab/doc/
  x = gpml_randn(0.8, 20, 1);                 % 20 training inputs
  y = sin(3*x) + 0.1*gpml_randn(0.9, 20, 1);  % 20 noisy training targets
  xs = linspace(-3, 3, 61)';                  % 61 test inputs 
end


function [x, y, xs] = gen_data2()
  % More point similar to first example from GPML tutorial:
  % http://www.gaussianprocess.org/gpml/code/matlab/doc/
  x = gpml_randn(0.8, 1000, 1);                 % training inputs
  y = sin(3*x) + 0.1*gpml_randn(0.9, 1000, 1);  % noisy training targets
  xs = linspace(-3, 3, 61)';                  % 61 test inputs 
end


function [locox] = locoTransformOne(x)
  if x < 0, locox = x.*20; else locox = x; end
end


function [] = plot_gp(mu, s2, x, y, xs, xu)
  % Simple plot of GP posterior.
  figure
  f = [mu+sqrt(s2); flipdim(mu-sqrt(s2),1)];
  fill([xs; flipdim(xs,1)], f, [255,234,253]/255)
  hold on; plot(xs, mu, '--m'); plot(x, y, '+b');
  if size(xu,1) > 0
      hold on; scatter(xu, ones(size(xu)), 'ok');
  end
end

