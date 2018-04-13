function swing_trajectory_generator(x0, dx0, x1, dx1, y0, dy0, y1, dy1, y_apex, T)

tf = T/2;
[a0 a1 a2 a3] = ThirdOrderTrajectory(y0, dy0, y_apex, 0.0, tf);
t = 0:0.001:tf;
ys1 = a3*t.^3 + a2*t.^2 + a1*t + a0;
yds1 = 3*a3*t.^2 + 2*a2*t + a1;

[a0 a1 a2 a3] = ThirdOrderTrajectory(y_apex, 0.0, y1, dy1, tf);
t = 0:0.001:tf;
ys2 = a3*t.^3 + a2*t.^2 + a1*t + a0;
yds2 = 3*a3*t.^2 + 2*a2*t + a1;

t = 0:0.001:T;
y_swing = [ys1 ys2(2:end)];
yd_swing = [yds1 yds1(2:end)];

[a0 a1 a2 a3] = ThirdOrderTrajectory(x0, dx0, x1, dx1, T);

x_swing = a3*t.^3 + a2*t.^2 + a1*t + a0;
xd_swing = 3*a3*t.^2 + 2*a2*t + a1;

swing_traj = [t' x_swing' xd_swing' y_swing' yd_swing'];
save('swing_traj.mat', 'swing_traj')
end