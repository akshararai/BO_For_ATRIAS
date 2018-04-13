% Six Hump Camel Back: 6 local minima, two of which are global.
% Global min: f(x,y) = -1.0316 at (-0.0898, 0.7126) and (0.0898, -0.7126)
% https://www.sfu.ca/~ssurjano/camel6.html
% https://github.com/SheffieldML/GPyOpt/tree/master/examples/six_hump_camel
function [res, out1, out2, out3, out4] = testf(X, in1, in2, in3)
x1 = X(1);
x2 = X(2);

term1 = (4-2.1*x1^2+(x1^4)/3) * x1^2;
term2 = x1*x2;
term3 = (-4+4*x2^2) * x2^2;

res = term1 + term2 + term3;

out1 = 0; out2 = 0; out3 = 0; out4 = 0;
end
