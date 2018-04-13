function [m,dm] = meanMmap(hyp, x)

% Similar to gpml_2017/mean/meanConst.m and meanZero.m
%
if nargin<2, m = '0'; return; end             % report number of hyperparameters
m = mmap(x);                                     % query the map to get the mean
dm = @(q) zeros(0,1);                                   % directional derivative
