function z = inverse_dynamics_double_support(GRFx1, GRFx2, GRFy1, GRFy2, theta, phi1, phi2, l_leg1, l_leg2, w_primary, w_secondary,...
                                  tau_rotor_f1, tau_rotor_f2, tau_rotor_b1, tau_rotor_b2, ...
                                  d_vertical_com, LEG_MTR_GEAR_RATIO, g,...
                                  j_leg, i_robot, ...
                                  m_total_real, m_leg, x_dx, y_dy)  
% syms g g_reduced N w real
% syms t_rf1 t_rb1 t_rf2 t_rb2 real
% syms m m_leg Il It d real
% syms t_sf1 t_sb1 t_sf2 t_sb2 real
% syms l1 g1 l2 g2 dl1 dg1 dl2 dg2 ddl1 ddg1 ddl2 ddg2 real
% syms GRFx1 GRFy1 GRFx2 GRFy2 real
% syms x y th dx dy dth ddx ddy ddth real
% 

%z =[ddx ddy ddl1 ddl2 ddg1 ddg2 ddth t11 t12 t01 t02];
%Parameters
    m      = m_total_real;    % mass of system + boom
    Il     = j_leg;             % inertia of leg
    It     = i_robot;           % inertia of system3.
    d      = d_vertical_com;             % distance from hip to system COM
    g1     = phi1(1)+theta(1);   % global leg angle (after spring)
    dg1    = phi1(2)+theta(2);
    g2     = phi2(1)+theta(1);
    dg2    = phi2(2)+theta(2);
    
    th     = theta(1);          % global torso angle   
    dth    = theta(2);
    t_rf1   = tau_rotor_f1;       % torque applied at front rotor 
    t_rb1   = tau_rotor_b1;       % torque applied at back rotor
    t_rf2   = tau_rotor_f2;       % torque applied at front rotor 
    t_rb2   = tau_rotor_b2;       % torque applied at back rotor

    N      = -LEG_MTR_GEAR_RATIO;
    l1      = l_leg1(1);          % leg length
    l2      = l_leg2(1);          % leg length
    dl1     = l_leg1(2);
    dl2     = l_leg2(2);

    dx = x_dx(2);
    dy = y_dy(2);
    
    J = [1 0 sin(g1) 0 l1*cos(g1) 0 -d*cos(th);
    0 1 cos(g1) 0 -l1*sin(g1) 0 d*sin(th);
    1 0 0 sin(g2) 0 l2*cos(g2) -d*cos(th);
    0 1 0 cos(g2) 0 -l2*sin(g2) d*sin(th)];

    dJ = [0 0 dg1*cos(g1) 0 -l1*dg1*sin(g1)+dl1*cos(g1) 0 d*dth*sin(th);
    0 0 -dg1*sin(g1) 0 -l1*dg1*cos(g1)-dl1*sin(g1) 0 d*dth*cos(th);
    0 0 0 dg2*cos(g2) 0 -l2*dg2*sin(g2)+dl2*cos(g2) d*dth*sin(th);
    0 0 0 -dg2*sin(g2) 0 -l2*dg1*cos(g1)-dl2*sin(g2) d*dth*cos(th)];

    F = [GRFx1 GRFy1 GRFx2 GRFy2]';
    dq = [dx dy dl1 dl2 dg1 dg2 dth]';

    M = [m 0 0 0 0 0 0;
    0 m 0 0 0 0 0;
    0 0 m_leg 0 0 0 0;
    0 0 0 m_leg 0 0 0;
    0 0 0 0 Il 0 0;
    0 0 0 0 0 Il 0;
    0 0 0 0 0 0 It];

    h = [0; -m*g; 0; 0; 0; 0; (-t_rf1-t_rf2-t_rb1-t_rb2)];

    S = [0 0 1/w_primary 0 1 0 -(N-1)/N;
    0 0 -1/w_primary 0 1 0 -(N-1)/N;
    0 0 0 1/w_secondary 0 1 -(N-1)/N;
    0 0 0 -1/w_secondary 0 1 -(N-1)/N];
    S = S';

    P = [M -S;
    J zeros(4,4)];

    H = [h;-dJ*dq];
    lambda = [J zeros(4,4)]'*F;

    z = P\(H+lambda);

end