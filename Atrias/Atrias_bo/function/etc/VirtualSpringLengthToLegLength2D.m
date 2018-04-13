%% Computes state of ATRIAS' primary leg given a desired virtual spring length and leg angle
% Only valid for 2D
% phi: Angle between torso and leg in 2D
% l_virtual: Distance between foot point and center of mass
% pelvis_to_com: Distance between pelvis and com
function [ l_leg, dl_leg ] = VirtualSpringLengthToLegLength2D( phi, dphi, l_virtual, dl_virtual, pelvis_to_com, d_pelvis_to_com )

    d = pelvis_to_com;
    dd = d_pelvis_to_com;
    
    l_leg = (2.*d.*cos(phi) + sqrt((2.*d.*cos(phi)).^2 + 4.*(l_virtual.^2 - d.^2))) ./ 2;
    dl_leg = dd.*cos(phi) - (2.*dd.*d - 2.*dl_virtual.*l_virtual - 2.*dd.*cos(phi).^2.*d + 2.*dphi.*cos(phi).*sin(phi).*d.^2)/(2.*(cos(phi).^2.*d.^2 - d.^2 + l_virtual.^2).^(1/2)) - dphi.*sin(phi).*d; 

end

