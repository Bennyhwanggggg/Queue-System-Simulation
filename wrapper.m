%% Wrapper script
clear
clc
ntests = importdata('num_tests.txt');

% rng settings
% Save
% rand_setting = rng;
% save rand_setting_4

% Load
% load rand_setting_2
% rng(rand_setting)

for i = 1:ntests
    display(i);
    mode = importdata(strcat('mode_', num2str(i), '.txt'));
    mode = mode{1};
    para = importdata(strcat('para_', num2str(i), '.txt'));
    m = para(1);
    setup_time = para(2);
    Tc = para(3);
    
    arrivals = importdata(strcat('arrival_', num2str(i), '.txt'));
    services = importdata(strcat('service_', num2str(i), '.txt'));
    %%%% Setup for Tc design task uncomment if required%%%%
%     mode = 'random';
%     m = 5;
%     setup_time = 5;
%     arrivals = 0.35;
%     services = 1;
%     Tc = 0.1;
    
    if strcmp(mode,'trace')
        time_end = 200;%-1;
        list_of_arrival_times = arrivals;
        list_of_service_times = services;
        [mrt, departure_times] = simulation(mode, list_of_arrival_times, list_of_service_times, m, setup_time, Tc, time_end);
    elseif strcmp(mode,'random')
        time_end = para(4);
%         time_end = 5000;
        lambda = arrivals;
        mu = services;
        [mrt, departure_times] = simulation(mode, lambda, mu, m, setup_time, Tc, time_end);
    end
    display(mrt)
%     display(departure_times)
    mrtfile = fopen(strcat('mrt_', num2str(i), '.txt'), 'w');
    fprintf(mrtfile, '%.3f', mrt);
    fclose(mrtfile);
    departurefile = fopen(strcat('departure_', num2str(i), '.txt'), 'w');
    for row = 1:size(departure_times,1)
        fprintf(departurefile,'%.3f\t %.3f\n', departure_times(row,:));
    end
    fclose(departurefile);
end