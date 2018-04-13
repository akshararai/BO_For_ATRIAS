%% ThirdOrderTrajectory    Solves coefficients for optimal third order trajectory polynomial with boundary constraints.
%
%    Lengths: meters
%    Angles:  radians
%    Time:    seconds
%    
%    [a0, a1, a2, a3] = ThirdOrderTrajectory(x0, dx0, xf, dxf, tf) determines the four coefficients for a third order
%    polynomial position trajectory that minimizes torque squared. The trajectory starts at x0 with dx0 velocity and
%    ends at xf with dxf velocity. The trajectory executes in tf seconds.
%
%    position(t) = a3*t^3 + a2*t^2 + a1*t + a0;
%    velocity(t) = 3*a3*t^2 + 2*a2*t + a1;
%    acceleration(t) = 6*a3*t + 2*a2;
function [a0, a1, a2, a3] = ThirdOrderTrajectory(x0, dx0, xf, dxf, tf)
    a0 = x0;
    a1 = dx0;
    a2 = -(3*x0 - 3*xf + 2*dx0*tf + dxf*tf)/tf^2;
    a3 = (2*x0 - 2*xf + dx0*tf + dxf*tf)/tf^3;
end