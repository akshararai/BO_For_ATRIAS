%% FifthOrderTrajectory    Solves coefficients for fifth order trajectory polynomial with boundary constraints.
%
%    Lengths: meters
%    Angles:  radians
%    Time:    seconds
%    
%    [a0, a1, a2, a3] = FifthOrderTrajectory(x0, dx0, ddx0, xf, dxf, ddxf, tf) determines the five coefficients for a
%    fifth order polynomial position trajectory. Note that this is not optimal for torque squared.
%    The trajectory starts at (x0, dx0, ddx0) and ends at (xf, dxf, ddxf). The trajectory executes in tf seconds.
%
%    position(t) = a5*t^5 + a4*t^4 + a3*t^3 + a2*t^2 + a1*t + a0;
%    velocity(t) = 5*a5*t^4 + 4*a4*t^3 + 3*a3*t^2 + 2*a2*t + a1;
%    acceleration(t) = 20*a5*t^3 + 12*a4*t^2 + 6*a3*t + 2*a2;
function [a0, a1, a2, a3, a4, a5] = FifthOrderTrajectory(x0, dx0, ddx0, xf, dxf, ddxf, tf)
    a0 = x0;
    a1 = dx0;
    a2 = ddx0/2;
    a3 = -(20*x0 - 20*xf + 12*dx0*tf + 8*dxf*tf + 3*ddx0*tf^2 - ddxf*tf^2)/(2*tf^3);
    a4 = (30*x0 - 30*xf + 16*dx0*tf + 14*dxf*tf + 3*ddx0*tf^2 - 2*ddxf*tf^2)/(2*tf^4);
    a5 = -(12*x0 - 12*xf + 6*dx0*tf + 6*dxf*tf + ddx0*tf^2 - ddxf*tf^2)/(2*tf^5);
end