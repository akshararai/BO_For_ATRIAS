% Plot AMTI force plate data from text file 
function [ F , t_force_plate] = PlotForcePlateData( file_name , time_offset, time_scale)
    
    if nargin < 2
        time_offset = 0;
        time_scale = 0.9982;
    end
    if nargin < 3
        time_scale = 0.9982;
    end

    force_plate_data = load(file_name);
    N_force_plate_data = size(force_plate_data,1);
    freq_force_plate = 200;
    t_force_plate = 0:(1/freq_force_plate):(N_force_plate_data/freq_force_plate)-(1/freq_force_plate);
    t_force_plate = (time_offset + t_force_plate*time_scale)';
    figure('Name','Force plate data');
    subplot(211);
    title('Vertical Force');
    plot(t_force_plate, force_plate_data(:,3));
    legend('F_z');
    subplot(212);
    title('Transverse Forces');
    plot(t_force_plate, force_plate_data(:,1:2));
    legend('F_x','F_y');

    F = force_plate_data(:,1:3);
end