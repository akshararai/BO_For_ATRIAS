%% Calculate corresponding leg angles and leg length for a 3D COM position
%
% Outputs:
% All outputs are of the form [position, velocity, acceleration]
%     Alpha:            relative leg sagittal angle 
%     Beta:             relative leg lateral angle
%     LegLength:        leg length
%
% Inputs:
%     X:                x position of COM relative to foot
%     Y:                y position of COM relative to foot
%     Z:                z position of COM relative to foot
%     right_left:       true if using right leg, false if using left leg
%     com_pelvis_body:  Pelvis to COM distance in body frame
%     dcom_pelvis_body: Pelvis to COM velocity in body frame
%     q:              position vector with the following order
%                         roll
%                         pitch
%                         relative leg sagittal angle 
%                         relative leg lateral angle 
%                         leg length

function [Alpha, Beta, LegLength]  =  ...
          ConvertComPositionToLegAnglesAndLegLength(X, Y, Z, right_left,...
                                                     com_pelvis_body, dcom_pelvis_body,...
                                                     q, dq, ...
                                                     lateral_offset)                                             
                                              
    phi         = q(1); % roll
    theta       = q(2); % pitch
    alpha       = q(3); % sagittal angle
    beta        = q(4); % lateral angle
    l_leg       = q(5); % leg length
    dphi        = dq(1); 
    dtheta      = dq(2); 
    dalpha      = dq(3); 
    dbeta       = dq(4); 
    dl_leg      = dq(5);

    % map this vector to the body-frame, also compute com-pelvis such that we
    % can solve for pelvis-foot in body frame, which then maps through IK to
    % joint angles and velocities.
    R_by = R_roll(phi)'*R_pitch(theta)'; %from yawed frame to body frame
    dR_by = dR_roll(phi,dphi)'*R_pitch(theta)' + R_roll(phi)'*dR_pitch(theta,dtheta)';
    
    % com-pelvis = pelvis2com in body frame (current state)
    ddp2c_b = [0;0;0]; % ignore acceleration of pelvis to com in body frame

    %% Left Leg
    if ~right_left
        % com-foot = foot2com in body frame (target value)
        fL2c_b = R_by*[X(1); Y(1); Z(1)]; 
        dfL2c_b = R_by*[X(2); Y(2); Z(2)] + dR_by*[X(1); Y(1); Z(1)]; 
        ddfL2c_b = R_by*[X(3); Y(3); Z(3)] + 2*dR_by*[X(2); Y(2); Z(2)];
        % pelvis-foot = foot2pelvis in body frame (target value)
        fL2p_b = fL2c_b - com_pelvis_body; % target position
        dfL2p_b = dfL2c_b - dcom_pelvis_body; % target velocity
        ddfL2p_b = ddfL2c_b - ddp2c_b; % target accel
        [alpha_star_L, beta_star_L, l_star_L] = left_leg_IK(-fL2p_b, lateral_offset); %IK wants foot relative to pelvis
        [dalpha_star_L, dbeta_star_L, dl_star_L] = velocity_map_L(alpha, beta, l_leg, ... %use current robot state
                               dfL2p_b(1), dfL2p_b(2), dfL2p_b(3), lateral_offset); % derived jacobian uses velocity of pelvis w.r.t foot.
        [ddalpha_star_L, ddbeta_star_L, ddl_star_L] = accel_map_L(alpha, beta, l_leg, dalpha, dbeta, dl_leg, ...
            ddfL2p_b(1), ddfL2p_b(2), ddfL2p_b(3), lateral_offset);                                          

        % unwrap alpha_star to [0 2*pi]
        while alpha_star_L < 0
          alpha_star_L = alpha_star_L + 2*pi;
        end
        while alpha_star_L > 2*pi
          alpha_star_L = alpha_star_L - 2*pi;
        end
        
        Alpha = [alpha_star_L; dalpha_star_L; ddalpha_star_L];
        Beta = [beta_star_L; dbeta_star_L; ddbeta_star_L];
        LegLength = [l_star_L; dl_star_L; ddl_star_L];  

    %% Right Leg
    else
        % com-foot = foot2com in body frame (target value)
        fR2c_b = R_by*[X(1); Y(1); Z(1)]; 
        dfR2c_b = R_by*[X(2); Y(2); Z(2)] + dR_by*[X(1); Y(1); Z(1)]; 
        ddfR2c_b = R_by*[X(3); Y(3); Z(3)] + dR_by*[X(2); Y(2); Z(2)];
        % pelvis-foot = foot2pelvis in body frame (target value)
        fR2p_b = fR2c_b - com_pelvis_body; % target position
        dfR2p_b = dfR2c_b - dcom_pelvis_body; % target velocity
        ddfR2p_b = ddfR2c_b - ddp2c_b; % target accel
        [alpha_star_R, beta_star_R, l_star_R] = right_leg_IK(-fR2p_b, lateral_offset); %IK wants foot relative to pelvis
        [dalpha_star_R, dbeta_star_R, dl_star_R] = velocity_map_R(alpha, beta, l_leg, ... %use current robot state
                               dfR2p_b(1), dfR2p_b(2), dfR2p_b(3), lateral_offset); % derived jacobian uses velocity of pelvis w.r.t foot.
        [ddalpha_star_R, ddbeta_star_R, ddl_star_R] = accel_map_R(alpha, beta, l_leg, dalpha, dbeta, dl_leg, ...
            ddfR2p_b(1), ddfR2p_b(2), ddfR2p_b(3), lateral_offset);                                          

        % unwrap alpha_star to [0 2*pi]
        while alpha_star_R < 0
          alpha_star_R = alpha_star_R + 2*pi;
        end
        while alpha_star_R > 2*pi
          alpha_star_R = alpha_star_R - 2*pi;
        end

        Alpha = [alpha_star_R; dalpha_star_R; ddalpha_star_R];
        Beta = [beta_star_R; dbeta_star_R; ddbeta_star_R];
        LegLength = [l_star_R; dl_star_R; ddl_star_R];  
    end

end

function R_yp = R_pitch(theta) % map from pitched to yawed
  R_yp = [cos(theta), 0, sin(theta); 0 1 0; -sin(theta), 0, cos(theta)];
end

function dR_yp = dR_pitch(theta,dtheta)
  dR_yp = dtheta*[-sin(theta), 0, cos(theta); 0 0 0; -cos(theta), 0, -sin(theta)];
end

function R_pr = R_roll(phi) % map from rolled to pitched
  R_pr = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)]; 
end

function dR_pr = dR_roll(phi,dphi)
  dR_pr = dphi*[0 0 0; 0 -sin(phi) -cos(phi); 0 cos(phi) -sin(phi)];
end

function [alpha1, beta1, l1] = left_leg_IK(pelvis2foot, r)
  xf = pelvis2foot(1); yf = pelvis2foot(2); zf = pelvis2foot(3);
  l1 = sqrt(max(0,xf^2+yf^2+zf^2-r^2));
  alpha1 = atan2(xf, -sqrt(max(0,l1^2-xf^2))); % assumes alpha1 > 90
  beta1_a = 2*atan((yf-sqrt(max(0,yf^2+zf^2-r^2)))/ (r+zf));
  beta1_b = 2*atan((yf+sqrt(max(0,yf^2+zf^2-r^2)))/ (r+zf));
  error_a = (r*sin(beta1_a)-l1*cos(alpha1)*cos(beta1_a) - yf)^2 + ...
            (r*cos(beta1_a)+l1*cos(alpha1)*sin(beta1_a) - zf)^2;
  error_b = (r*sin(beta1_b)-l1*cos(alpha1)*cos(beta1_b) - yf)^2 + ...
            (r*cos(beta1_b)+l1*cos(alpha1)*sin(beta1_b) - zf)^2;
  if error_a <= error_b
    beta1 = beta1_a;
  else
    beta1 = beta1_b;
  end
  if alpha1< 0
    alpha1 = alpha1+2*pi;
  end
end

% map between xyz velocities and abl velocities dx = J dalpha
function J = leftleg_J(alpha_L, beta_L, l_L, r)
  J = [ -l_L*cos(alpha_L),                                          0,            -sin(alpha_L);
        -l_L*cos(beta_L)*sin(alpha_L), - r*cos(beta_L) - l_L*cos(alpha_L)*sin(beta_L),  cos(alpha_L)*cos(beta_L);
         l_L*sin(alpha_L)*sin(beta_L),   r*sin(beta_L) - l_L*cos(alpha_L)*cos(beta_L), -cos(alpha_L)*sin(beta_L)];
end

function dJ = leftleg_dJ(alpha1, beta1, l1, dalpha1, dbeta1, dl1, r)
  dJ = [dalpha1*l1*sin(alpha1) - dl1*cos(alpha1), 0, -dalpha1*cos(alpha1);
        dbeta1*l1*sin(alpha1)*sin(beta1) - dalpha1*l1*cos(alpha1)*cos(beta1) - dl1*cos(beta1)*sin(alpha1), dbeta1*(r*sin(beta1) - l1*cos(alpha1)*cos(beta1)) - dl1*cos(alpha1)*sin(beta1) + dalpha1*l1*sin(alpha1)*sin(beta1), - dalpha1*cos(beta1)*sin(alpha1) - dbeta1*cos(alpha1)*sin(beta1);
        dl1*sin(alpha1)*sin(beta1) + dalpha1*l1*cos(alpha1)*sin(beta1) + dbeta1*l1*cos(beta1)*sin(alpha1), dbeta1*(r*cos(beta1) + l1*cos(alpha1)*sin(beta1)) - dl1*cos(alpha1)*cos(beta1) + dalpha1*l1*cos(beta1)*sin(alpha1),   dalpha1*sin(alpha1)*sin(beta1) - dbeta1*cos(alpha1)*cos(beta1)];
end

% conversion from body-frame velocities of pelvis-foot to joint angles
function [dalpha1, dbeta1, dl1] = velocity_map_L(alpha1, beta1, l1, dx, dy, dz, r)
  % [dx dy dz]_body = J_left * [dalpha1, dbeta1, dl];
  dvec = leftleg_J(alpha1,beta1,l1,r)\[dx;dy;dz];
  dalpha1 = dvec(1); dbeta1 = dvec(2); dl1 = dvec(3);
end

function [ddalpha1, ddbeta1, ddl1] = accel_map_L(alpha1,beta1,l1, dalpha1, dbeta1, dl1, ...
    ddx_star, ddy_star, ddz_star, r)
  ddvec = leftleg_J(alpha1,beta1,l1,r)\ ([ddx_star;ddy_star;ddz_star] ...
            - leftleg_dJ(alpha1,beta1,l1,dalpha1,dbeta1,dl1,r)*[dalpha1;dbeta1;dl1]);
  ddalpha1 = ddvec(1); ddbeta1 = ddvec(2); ddl1 = ddvec(3);
end

function [alpha1, beta1, l1] = right_leg_IK(pelvis2foot, r)
  [alpha1, beta1, l1] = left_leg_IK(diag([1 -1 1])*pelvis2foot, r);
end

function J = rightleg_J(alpha_R, beta_R, l_R, r)
  J = diag([1 -1 1])*leftleg_J(alpha_R, beta_R, l_R, r);
end

function dJ = rightleg_dJ(alpha1, beta1, l1, dalpha1, dbeta1, dl1, r)
  dJ = diag([1 -1 1])*leftleg_dJ(alpha1, beta1, l1, dalpha1, dbeta1, dl1, r);
end

function [dalpha1, dbeta1, dl1] = velocity_map_R(alpha1, beta1, l1, dx, dy, dz, r)
  dvec = rightleg_J(alpha1,beta1,l1,r)\[dx;dy;dz];
  dalpha1 = dvec(1); dbeta1 = dvec(2); dl1 = dvec(3);
end

function [ddalpha1, ddbeta1, ddl1] = accel_map_R(alpha1,beta1,l1, dalpha1, dbeta1, dl1, ...
    ddx_star, ddy_star, ddz_star, r)
  ddvec = rightleg_J(alpha1,beta1,l1,r)\ ([ddx_star;ddy_star;ddz_star] ...
            -rightleg_dJ(alpha1,beta1,l1,dalpha1,dbeta1,dl1,r)*[dalpha1;dbeta1;dl1]);
  ddalpha1 = ddvec(1); ddbeta1 = ddvec(2); ddl1 = ddvec(3);
end