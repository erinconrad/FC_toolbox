function [is_seeg_time,seeg_stage] = match_seeg_times(times,file_index,seeg_secs,stage)

is_seeg_time = zeros(size(times));
seeg_stage = cell(size(times));

for t = 1:length(times)

    if isnan(times(t,1)), continue; end

    if file_index(t) > 1
        break
    end

    curr_time = times(t);
    run_times = [curr_time-30,curr_time+30]; % time is middle of 60 second window

     % Figure out if first run time falls between two transitions
    run_start_minus_seeg = run_times(1)-seeg_secs;

    if all(run_start_minus_seeg>0) % this block is too late, stop loop
        break
    end

    if all(run_start_minus_seeg<0) % this block too early, go to next loop
        continue
    end

    % find the first seeg transition that comes BEFORE
    prior_transition = find(run_start_minus_seeg<0);
    prior_transition = prior_transition(1)-1;

    % confirm that run time end comes before the next one
    if run_times(2) > seeg_secs(prior_transition+1)
        continue
    end

    % If made it here, then the run time falls between two state
    % transitions, and so I should be able to compare
    is_seeg_time(t) = 1;
    seeg_stage(t) = stage(prior_transition);
end

is_seeg_time = logical(is_seeg_time);

end