%%
% COMP9334 Project
% Author: Kuan-Chun Hwang
%
% Inputs: mode: random/trace
%         arrival:
%               random: lambda -> The inter-arrival probability distribution is exponentially distributed with parameter . This
% means the mean arrival rate of the jobs is . You will need to supply the value of to your
% program using the input parameter arrival.
%               trace: list of arrival time
%
%         service:
%               random: mu -> Let sk denote the service time of the k-th job arriving at the dispatcher. Each sk is the sum
% of three random numbers s1k, s2k and s3k, i.e sk = s1k + s2k + s3k where s1k,
% s2k and s3k are exponentially distributed random numbers with parameter . You will need
% to supply the value of  to your program using the input parameter service.
%               trace: list of service time
%
%         m: number of servers. Positive int
%
%         setup_time: set up time for computers. Positive float
%
%         delayedoff_time: Init value of countdown timer Tc. Postive float
%
%         time_end: Time when master clock end. Only relevent when mode is
%         random.
%%
function [mrt, departure_times] = simulation(mode, arrival, service, m, setup_time, delayedoff_time, time_end)
lambda = arrival;
mu = service;

% Option to see arrival and service time distribution for random mode
displaydistributions = 0;

response_time_cumulative = 0; %  The cumulative response time
num_customer_served = 0; % number of completed customers at the end of the simulation
departure_times = [];
% Event arrays setup
next_departure_time = Inf * ones(m,1);
next_setup_time = Inf * ones(m,1);
next_delayed_off_time = Inf * ones(m,1);

% To store random variable result to check distribution
random_arrival_result = [];
service_result = [];

% server_status: 0=OFF, 1=SETUP, 2=DELAYEDOFF, 3=BUSY
master_clock = 0;
server_status = zeros(m,1);
arrival_time_next_departure = zeros(m,1); % job finish time

% Initialise buffer
buffer_content = [];
queue_length = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialising the events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialising the arrival event if random mode
if strcmp(mode, 'random')
    next_arrival_time = -log(1-rand(1))/lambda;
    service_time_next_arrival = -log(1-rand(1))/mu;
    condition = master_clock <= time_end;
elseif strcmp(mode, 'trace')
    arrival_n = 1;
    next_arrival_time = arrival(arrival_n);
    service_time_next_arrival = service(arrival_n);
    time_end = length(arrival);
    arrival_n = arrival_n + 1;
    condition = num_customer_served < time_end;
end

e1 = 1;
s1 = 1;

while condition
    % Find the server with the next event
    [first_departure_time,first_departure_server] = min(next_departure_time);
    [first_setup_time,first_setup_server] = min(next_setup_time);
    [first_delayoff_time,first_delayoff_server] = min(next_delayed_off_time);
    
    % Find out whether the next event is an arrival or depature
    %
    % Need to change here to account for the each server status change
    % server_status: 0=OFF, 1=BUSY, 2=SETUP, 3=DELAYEDOFF
    % get the min of each status array and change next event tupe
    % accordingly
    % Look at event time arrays
    if (next_arrival_time < first_departure_time)
        next_event_time = next_arrival_time;
        next_event_type = 1;
    else
        next_event_time = first_departure_time;
        % first departure server has already been found just now
        next_event_type = 0;
    end
    
    if (next_arrival_time < first_departure_time && next_arrival_time < first_setup_time && next_arrival_time < first_delayoff_time)
       next_event_time = next_arrival_time;
       next_event_type = 1;
    elseif (first_departure_time < next_arrival_time && first_departure_time < first_setup_time && first_departure_time < first_delayoff_time)
        next_event_time = first_departure_time;
        % first departure server has already been found just now
        next_event_type = 0;
    elseif (first_setup_time < next_arrival_time && first_setup_time < first_departure_time && first_setup_time < first_delayoff_time)
        next_event_time = first_setup_time;
        % first setup server has already been found just now
        next_event_type = 2;
    elseif (first_delayoff_time < next_arrival_time && first_delayoff_time < first_departure_time && first_delayoff_time < first_setup_time)
        next_event_time = first_delayoff_time;
        % first delayoff server has already been found just now
        next_event_type = 3;
    end
    % update master clock
    master_clock = next_event_time;

    % take actions depending on the event type
    if (next_event_type == 1) % an arrival
        % If no server in delayed_off, we check if any are currently off
        % and turn it on and add to dispatcher as marked, otherwise add it to dispatcher as unmarked
        if ~ismember(2, server_status)
            if (ismember(0, server_status))
                buffer_content = [buffer_content ; next_arrival_time service_time_next_arrival 1];
                idle_server = min(find(server_status == 0));
                server_status(idle_server) = 1;
                next_setup_time(idle_server) = master_clock+setup_time;
            else
                buffer_content = [buffer_content ; next_arrival_time service_time_next_arrival 0];
            end
            queue_length = queue_length + 1;
        % If at least one server is in delayed off, find the one that is
        % furtherest away from timing out (smallest time as it only just
        % started counting down) and change it back to busy
        else 
            % We need to find which server has the highest delayed off time
            if ismember(Inf, next_delayed_off_time)
                [largest_delayed_off_time, ~] = max(next_delayed_off_time(next_delayed_off_time<max(next_delayed_off_time)));
                largest_delayed_off_server = find(next_delayed_off_time == largest_delayed_off_time);
            else
                [~, largest_delayed_off_server] = max(next_delayed_off_time);
            end
            server_status(largest_delayed_off_server) = 3; % set to busy
            next_delayed_off_time(largest_delayed_off_server) = Inf; % server no longer in delayed off so change timer back to inf
            next_departure_time(largest_delayed_off_server) = next_arrival_time + service_time_next_arrival;
            % arrival time of next departure tells us the arrival time of
            % the job so we can use it for response time calculation after
            arrival_time_next_departure(largest_delayed_off_server) = next_arrival_time;
        end
        % Get next arrival time and service time
        if (strcmp(mode, 'random'))
            random_arrival = log(1-rand(1))/lambda;
            random_arrival_result(e1) = -random_arrival;
            e1 = e1 + 1;
            % next_arrival_time = master_clock - log(1-rand(1))/lambda;
            next_arrival_time = master_clock - random_arrival;
            service_time_next_arrival = (-log(1-rand(1))/mu)+(-log(1-rand(1))/mu)+(-log(1-rand(1))/mu);
            service_result(s1) = service_time_next_arrival;
            s1 = s1+1;
        elseif (strcmp(mode, 'trace'))
            if arrival_n <= length(arrival)
                next_arrival_time = arrival(arrival_n);
                service_time_next_arrival = service(arrival_n);
                arrival_n = arrival_n + 1;
            else
                next_arrival_time = Inf;
            end
        end
    elseif (next_event_type == 2) % a server has finished setup and is ready to take a job off queue
        if (queue_length && buffer_content(1, 3)==1)% If there's any job in the queue
            % Change the status of first setup server to busy
            server_status(first_setup_server) = 3;
            next_setup_time(first_setup_server) = Inf;
            % Schedule the next departure event using the first element in
            % the buffer (first row)
            next_departure_time(first_setup_server) = master_clock + buffer_content(1,2);
            arrival_time_next_departure(first_setup_server) = ...
                buffer_content(1,1);
            % remove customer from buffer and decrement queue length
            buffer_content(1,:) = [];
            queue_length = queue_length - 1;
        end
    elseif (next_event_type == 3) % A delayed off server has expired count down, the server goes off
        if (first_delayoff_time >= delayedoff_time)
            server_status(first_delayoff_server) = 0;
            next_delayed_off_time(first_delayoff_server) = Inf;
        end
    elseif (next_event_type == 0) % A job has finished in a server and will now depart
        response_time_cumulative = response_time_cumulative + master_clock - arrival_time_next_departure(first_departure_server);
        num_customer_served = num_customer_served + 1;
        departure_times = [departure_times; round(arrival_time_next_departure(first_departure_server), 3) round(master_clock, 3)];
        % If there are marked jobs waiting in queue, take the job,
        % otherwise go to delayed off
        if queue_length
            next_departure_time(first_departure_server) = ...
                master_clock + buffer_content(1,2);
            arrival_time_next_departure(first_departure_server) = ...
                buffer_content(1,1);
            buffer_content(1,:) = [];
            queue_length = queue_length - 1;
        else
            next_departure_time(first_departure_server) = Inf;
            server_status(first_departure_server) = 2;
            next_delayed_off_time(first_departure_server) = master_clock + delayedoff_time;
        end
    end
    
    number_of_setup_servers = sum(server_status(:) == 1);
    number_of_marked_jobs = sum(buffer_content(:, 3)==1);
    number_of_unmarked_jobs = sum(buffer_content(:, 3)==0);
    num_diff = number_of_setup_servers - number_of_marked_jobs;
    if num_diff > 0
        if number_of_unmarked_jobs < num_diff
            unmarked_job_rows = find(buffer_content(:, 3) == 0);
            % mark number of jobs equal amount of setup server
            for i=1:number_of_unmarked_jobs
                buffer_content(unmarked_job_rows(i), 3) = 1;
            end
            % turnoff additional server starting with the longest setup
            % time one
            number_of_server_to_turnoff = num_diff - number_of_unmarked_jobs;
            for i=1:number_of_server_to_turnoff
                if ismember(Inf, next_setup_time)
                    [val, ~] = max(next_setup_time(next_setup_time<max(next_setup_time)));
                    server_to_turn_off = find(next_setup_time == val);
                else
                    [~, server_to_turn_off] = max(next_setup_time);
                end
                if server_status(server_to_turn_off) == 1
                    next_setup_time(server_to_turn_off) = Inf;
                    server_status(server_to_turn_off)= 0;
                end
            end
        elseif number_of_unmarked_jobs >= num_diff % If there are more unmarked jobs then mark as many job as possible to marked
            unmarked_job_rows = find(buffer_content(:, 3) == 0);
            for i=1:num_diff  
                buffer_content(unmarked_job_rows(i), 3) = 1;
            end
        end
    end
    
    if strcmp(mode, 'random')
        condition = master_clock <= time_end;
    elseif strcmp(mode, 'trace')
        condition = num_customer_served < time_end;
    end
end


mrt = round(response_time_cumulative/num_customer_served, 3);

% If distribution graph required
if strcmp(mode, 'random') && displaydistributions
    figure(1)
    histogram(random_arrival_result, 50)
    title('Distribution of arrival times')
    xlabel('next arrival time')
    ylabel('frequency')
    figure(2)
    histogram(service_result, 50)
    title('Distribution of service times')
    xlabel('next service time')
    ylabel('frequency')
end
end











