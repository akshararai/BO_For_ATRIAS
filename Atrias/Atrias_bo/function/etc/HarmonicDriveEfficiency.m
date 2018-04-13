% Compute efficiency of harmonic drive at given speed (rpm) and load torque (Nm)
function [ mu ] = HarmonicDriveEfficiency( rpm, load_torque )
    % rpm = [500, 1000, 2000, 3500], mu = [0.7274, 0.68493, 0.62329, 0.56027]
    harmonic_efficiency_at_rated_torque = 9.422e-9*min(abs(rpm),3500).^2 - 9.296e-5*min(abs(rpm),3500) + 0.7704;
    % tau = [0.1:0.1:1], compensation = [0.3, 0.48657, 0.6194, 0.71791, 0.79254, 0.85075, 0.89851, 0.93881, 0.97313, 1]
    harmonic_efficiency_compensation_coefficient = (1.362*max(min(abs(load_torque/172),1),0.3)+0.00256)./(0.3636+max(min(abs(load_torque/172),1),0.3));
    mu = harmonic_efficiency_at_rated_torque.*harmonic_efficiency_compensation_coefficient;
end