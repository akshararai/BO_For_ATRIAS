% Calculate 4 spring SEA spring forces based on 4 deflections
% RB RF LB LF
function [ tau_sea ] = SeaSpringForce( deflections )

    direction = sign(deflections);
    x = min(abs(deflections),0.42);
    tau_sea = [-1960.6*x(1)^2 + 3579.2*x(1); ...
             -4954.9*x(2)^2 + 4334*x(2); ...
             -2836.1*x(3)^2 + 3698*x(3); ...
             -1846.9*x(4)^2 + 4044.4*x(4)];
    tau_sea = direction .* tau_sea;

end
