% Braning function.
% Global min: f(x,y) = 0.397887 
% x* =  (-pi, 12.275), (pi, 2.275), (9.42478, 2.475)
% x1 in [-5, 10], x2 in [0, 15]. 
% https://www.sfu.ca/~ssurjano/branin.html
function [res, out1, out2, out3, out4] = braninf(X, in1, in2, in3)
x1 = X(1);
x2 = X(2);

t = 1 / (8*pi);
s = 10;
r = 6;
c = 5/pi;
b = 5.1 / (4*pi^2);
a = 1;

term1 = a * (x2 - b*x1^2 + c*x1 - r)^2;
term2 = s*(1-t)*cos(x1);

res = term1 + term2 + s;
res = res + 1e-2*rand;

out1 = 0; out2 = 0; out3 = 0; out4 = 0;
end
