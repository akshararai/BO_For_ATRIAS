close all 
clear all

g = 9.81;
z = 0.8;

T = 1.0;

dt = 0.001;
time = 0:dt:T;

m_total_real = 62.5910;
i_robot = 4.133;
desired_torso_pitch = 0;
desired_z = 0.8;

A =  [0 0 1 0; % x th dx dth
      0 0 0 1; 
      0 0 0 0; 
      -m_total_real*g/i_robot 0 0 0];
B = [0; 0; 1/m_total_real; desired_z/i_robot];

A = A*dt + eye(size(A));
B = B*dt;

Q = diag([100 100 0.01 100]); 
R = 0.01;

x0 = [0.1; 20*pi/180; 0.0; 0.0];
x_des = [0.0; 0.0; 0.0; 0.0];
u_des = 0; %m_total_real*g; 

K = dlqr(A, B, Q, R);

x(1,:) = x0;
u(1) = -K*x(1,:)' + u_des;
for i = 1: length(time)-1
    x(i+1,:) = A*x(i,:)' + B*u(i); %    xd(i,:) = A*x(i,:)' + B*u(i); %x(i+1,:) = xd(i,:)*dt + x(i,:);
    u(i+1) = -K*(x(i+1,:)-x_des')' + u_des;
end
% 
figure(1), clf, hold on
title('COM position');
plot(x(:,1))

figure(2), clf, hold on
title('theta');
plot(x(:,2))

figure(3), clf, hold on
title('Fx');
plot(u)

err_lqr= norm(x)

%% FH LQR on balancing

H = diag([10000 10000 100 100]);
P_final = H;
K_final = (R + B'*P_final*B)\(B'*P_final*A);
P = zeros( 4,4,length(time));
K = zeros( 4, length(time));

P(:,:, end) = P_final;
K(:, end) = K_final;

for i = length(time)-1:-1:1
    P(:,:,i) = A'*P(:,:,i+1)*A - (A'*P(:,:, i+1)*B)*inv(R+B'*P(:,:,i+1)*B)*(B'*P(:,:,i+1)*A)+Q;
    K(:,i) = (R+B'*P(:,:,i)*B)\(B'*P(:,:,i)*A);
end
x(1,:) = x0;

for i = 1: length(time)-1
    x(i+1,:) = A*x(i,:)' + B*u(i);
    u(i+1) = -K(:,i)'*x(i+1,:)';
end

figure(4), clf, hold on
title('COM position');
plot(x(:,1))

figure(5), clf, hold on
title('theta');
plot(x(:,2))

figure(6), clf, hold on
title('Fx');
plot(u)

err_lqr_fh= norm(x)

%% LQR with plan
clear x
plan = [1.0 0.0
1.0 0.1
1.0 0.2];

T = sum(plan(:,1));
p_xs = plan(:,2);

g = 9.81;

dt = 0.001;
time = 0:dt:T;

m_total_real = 62.5910;
i_robot = 4.133;
desired_torso_pitch = 0;
desired_z = 0.8;

A =  [0 0 1 0; % x th dx dth
      0 0 0 1; 
      0 0 0 0; 
      -m_total_real*g/i_robot 0 0 0];
B = [0; 0; 1/m_total_real; desired_z/i_robot];

A = A*dt + eye(size(A));
B = B*dt;

X_star = zeros(length(time), 4);

time_last = 0;
step_last = 1;
for i =1 :length(time)
    if(time(i)-time_last <= plan(step_last,1))
        X_star(i,1) = p_xs(step_last);
    else
        time_last = time(i);
        step_last = step_last + 1;
        X_star(i,1) = p_xs(step_last);
    end
end

Q = diag([100000 5 10 1]); 
R = 0.001;

x0 = [0.0; 0*pi/180; 0.0; 0.0];
x_des = [0.0; 0.0; 0.0; 0.0];
u_des = 0; %m_total_real*g; 


K = dlqr(A, B, Q, R);
x(1,:) = x0;
u(1) = -K*x(1,:)' + u_des;

for i = 1: length(time)-1
    x(i+1,:) = A*x(i,:)' + B*u(i);
    u(i+1) = -K*(x(i+1,:)'-X_star(i+1,:)');
end

K_lqr = K;
figure(1), clf, hold on
title('COM position lqr');
plot(x(:,1))
plot(X_star(:,1), 'r')

figure(2), clf, hold on
title('theta lqr');
plot(x(:,2))

figure(3), clf, hold on
title('Fx lqr');
plot(u)

err_lqr= norm(x-X_star)

%% Finite horizon

H = Q;
P_final = H;
K_final = (R + B'*P_final*B)\(B'*P_final*A);
P = zeros( 4,4,length(time));
K = zeros( 4, length(time));

P(:,:, end) = P_final;
K(:, end) = K_final;

for i = length(time)-1:-1:1
    P(:,:,i) = A'*P(:,:,i+1)*A - (A'*P(:,:, i+1)*B)*inv(R+B'*P(:,:,i+1)*B)*(B'*P(:,:,i+1)*A)+Q;
    K(:,i) = (R+B'*P(:,:,i)*B)\(B'*P(:,:,i)*A);
end

for i = 1: length(time)-1
    x(i+1,:) = A*x(i,:)' + B*u(i);
    u(i+1) = -K(:,i)'*(x(i+1,:)'-X_star(i+1,:)');
end

figure(4), clf, hold on
title('COM position fhlqr');
plot(x(:,1))
plot(X_star(:,1), 'r')

figure(6), clf, hold on
title('theta fhlqr');
plot(x(:,2))

figure(5), clf, hold on
title('Fx fhlqr');
plot(u)

err_lqr_fh= norm(x-X_star)

%% DDP
K = K_lqr;
U_star = zeros(length(time));

for iter = 1:100000000
    V_x(length(time),:) = 2*(x(end,:)-X_star(end,:))*Q;
    V_xx(:,:,length(time)) = 2*Q;
    
    delU(1) = 0;

    for i=length(time)-1:-1:1
        Z_x = V_x(i+1,:)*A + (x(i+1,:)-X_star(i+1,:))*Q;
        Z_u = V_x(i+1,:)*B + (u(i+1)-U_star(i+1))*R;
        Z_xx = A'*V_xx(:,:,i+1)*A + Q;
        Z_ux = B'*V_xx(:,:,i+1)*A;
        Z_uu = B'*V_xx(:,:,i+1)*B + R;
        K(i+1,:) = Z_uu\Z_ux;
        delU(i+1) = (Z_uu\Z_u')';
        V_x(i,:) = Z_x - Z_u * K(i+1,:);
        V_xx(:,:,i)= Z_xx - Z_ux'*K(i+1,:);
    end

    xnew=x;
    unew=u;
    for i=1:length(time)-1
        xnew(i+1,:) = (A*xnew(i,:)' + B*unew(i));
        unew(i+1) = -K(i+1,:)*(xnew(i+1,:)'-x(i+1,:)')  + u(i+1) - 0.5*delU(i+1);%(u(i+1,:)-0.05*delU(i+1,:)) - (K(:,:,i+1)*(xnew(i+1,:)'-x(i+1)'))';%-(K(:,:,i+1)*(xnew(i+1,:)'-X_star(i+1,:)'))'+ U_star(i+1,:);%+ U_star(i+1,:);%-(K(:,:,i+1)*(xnew(i+1,:)'-x(i+1,:)'))';
    end


    err_new= norm(xnew-X_star);
    err_old=norm(x-X_star);

    if norm(err_new - err_old ) < 10^-6
        break
    end
    x=xnew;
    u=unew;
    % pause
end
figure(7), clf, hold on
title('COM position ddp');
plot(x(:,1))
plot(X_star(:,1), 'r')

figure(8), clf, hold on
title('theta ddp');
plot(x(:,2))
plot(X_star(:,2), 'r')

figure(9), clf, hold on
title('Fx ddp');
plot(u)


err_ddp = err_new
com_traj_x = [time' x(:,1)];
save('com_traj_x.mat', 'com_traj_x')

% %% Full LQR
% clear x Q R
% 
% plan = [1.0 0.0
% 1.0 0.1
% 1.0 0.2];
% 
% T = sum(plan(:,1));
% p_xs = plan(:,2);
% 
% g = 9.81;
% 
% dt = 0.001;
% time = 0:dt:T;
% 
% m_total_real = 62.5910;
% i_robot = 4.133;
% desired_torso_pitch = 0;
% desired_z = 0.8;
% 
% Fx0 = 0.0;
% Fy0 = m_total_real*g;
% y0 = desired_z;
% Q = diag([100 1 100 1 100 1]);
% R = diag([0.001 0.001]);
% 
% H = Q;
% A = [0 1 0 0 0 0
%      0 0 0 0 0 0
%      0 0 0 1 0 0
%      0 0 0 0 0 0
%      0 0 0 0 0 1
%      -Fy0/I 0 Fx0/I 0 0 0];
% B = [0; m; 0; m; 0; -y0/I + x0/I];
% 
% A = A*dt + eye(size(A));
% B = B*dt;
% 
% P_final = H;
% K_final = (R + B'*P_final*B)\(B'*P_final*A);
% P = zeros( 6, 6,length(time));
% K = zeros( 6, length(time));
% x = zeros( length(time), 6);
% u = zeros( length(time), 1);
% 
% 
% P(:,:, end) = P_final;
% K(:, end) = K_final;
% 
% for i = length(time)-1:-1:1

% 
%     P(:,:,i) = A'*P(:,:,i+1)*A - (A'*P(:,:, i+1)*B)*inv(R+B'*P(:,:,i+1)*B)*(B'*P(:,:,i+1)*A)+Q;
%     K(:,i) = (R+B'*P(:,:,i)*B)\(B'*P(:,:,i)*A);
% end
% x(1,:) = x0;
% u(1) = -K(:,1)'*(x(1,:)' - x_des) + u_des;
% for i = 1: length(time)-3
%     A =  [0 0 1 0; % x th dx dth
%           0 0 0 1; 
%           0 0 0 0; 
%           -Fx_traj(i)/i_robot 0 0 0];
%     B = [0; 0; 1/m_total_real; x_traj(i)/i_robot];
% 
%     A = A*dt + eye(size(A));
%     B = B*dt;
% 
%     x(i+1,:) = A*x(i,:)' + B*u(i);
%     u(i+1) = -K(:,i+1)'*(x(i+1,:)' - x_des) + u_des;
% end
% 
% figure(3), clf, hold on
% title('COM x position')
% plot(x_traj)
% 
% figure(4), clf, hold on
% title('COM z position');
% plot(x(:,1))
% 
% figure(5), clf, hold on
% title('theta');
% plot(x(:,2))
% 
% figure(6), clf, hold on
% title('Fz');
% plot(u)


 

%% LQR over z-theta
% 
% %x_plan
% 
% plan = [1.0 0.0
% 1.0 0.2
% 1.0 0.4];
% 
% load com_traj_x.mat
% 
% x_traj = com_traj_x(:,2);
% xdd_traj = diff(diff(x_traj));
% 
% Fx_traj = m_total_real*xdd_traj;
% 
% T = sum(plan(:,1));
% dt = 0.001;
% time = 0:dt:T;
% 
% 
% x0 = [0.7; 10*pi/180; 0.0; 0.0];
% x_des = [0.8; 0; 0; 0];
% u_des = [m_total_real*g];
% Q = diag([1000000 10000 100 100]);  
% R = 0.001;
% 
% H = Q;
% A =  [0 0 1 0; % z th dz dth
%       0 0 0 1; 
%       0 0 0 0; 
%       -Fx_traj(end)/i_robot 0 0 0];
% B = [0; 0; 1/m_total_real; x_traj(end)/i_robot];
% 
% A = A*dt + eye(size(A));
% B = B*dt;
% 
% P_final = H;
% K_final = (R + B'*P_final*B)\(B'*P_final*A);
% P = zeros( 4,4,length(time)-2);
% K = zeros( 4, length(time)-2);
% x = zeros( length(time)-2, 4);
% u = zeros( length(time)-2, 1);
% 
% 
% P(:,:, end) = P_final;
% K(:, end) = K_final;
% 
% for i = length(time)-3:-1:1
%     A =  [0 0 1 0; % x th dx dth
%           0 0 0 1; 
%           0 0 0 0; 
%           -Fx_traj(i)/i_robot 0 0 0];
%     B = [0; 0; 1/m_total_real; x_traj(i)/i_robot];
% 
%     A = A*dt + eye(size(A));
%     B = B*dt;
% 
%     P(:,:,i) = A'*P(:,:,i+1)*A - (A'*P(:,:, i+1)*B)*inv(R+B'*P(:,:,i+1)*B)*(B'*P(:,:,i+1)*A)+Q;
%     K(:,i) = (R+B'*P(:,:,i)*B)\(B'*P(:,:,i)*A);
% end
% x(1,:) = x0;
% u(1) = -K(:,1)'*(x(1,:)' - x_des) + u_des;
% for i = 1: length(time)-3
%     A =  [0 0 1 0; % x th dx dth
%           0 0 0 1; 
%           0 0 0 0; 
%           -Fx_traj(i)/i_robot 0 0 0];
%     B = [0; 0; 1/m_total_real; x_traj(i)/i_robot];
% 
%     A = A*dt + eye(size(A));
%     B = B*dt;
% 
%     x(i+1,:) = A*x(i,:)' + B*u(i);
%     u(i+1) = -K(:,i+1)'*(x(i+1,:)' - x_des) + u_des;
% end
% 
% figure(3), clf, hold on
% title('COM x position')
% plot(x_traj)
% 
% figure(4), clf, hold on
% title('COM z position');
% plot(x(:,1))
% 
% figure(5), clf, hold on
% title('theta');
% plot(x(:,2))
% 
% figure(6), clf, hold on
% title('Fz');
% plot(u)
% 
% z_and_theta = [time(1:end-2)' x(:,1) x(:,2)];
% save('z_and_theta.mat','z_and_theta')
% err_lqr_fh= norm(x)
