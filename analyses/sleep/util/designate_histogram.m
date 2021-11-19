function [transitions,bins] = designate_histogram(sleep,n_periods,min_same,...
    later_search,time_to_take_spikes,rm_cluster)

% the transition is the first sleep period

%{
Take a list of binary sleep or no sleep labels. Loop over the indices. Look
back 3 hours and forward 3 hours to see if there are enough that are mostly
awake before and mostly asleep after. If it meets this criterion, it is a
potential sleep transition point. I then look 2 hours forward to see if
there is a better sleep transition point (one with more sleep where there
should be sleep and wake where there should be wake). I take the best one.
I then continue looping forward to build up a bunch of transition points.
%}

transitions = [];
i = time_to_take_spikes + 1;

while i < length(sleep) - time_to_take_spikes
    before = i-n_periods:i-1;
    after = i:i+n_periods;
    
    % if enough before are awake and 
    if sum(sleep(before) == 0)/n_periods > min_same && ...
            sum(sleep(after) == 1)/n_periods > min_same
        
        sum_appropriate = zeros(later_search+1,1);
        % this is a potential transition point, search later to see if
        % there is a better one
        for later = 0:later_search % up to 120 minutes later
            
            % count matching sleep/wake for this transition point
            new_before = i+later - n_periods:i+later-1;
            new_after = i+later:i+later+n_periods;
            
            sum_appropriate(later+1) = sum(sleep(new_before) == 0) + ...
                sum(sleep(new_after) == 1);
            
        end
        
        % find the max sum appropriate
        [n_match,I] = max(sum_appropriate);
        
        true_later = I;
        transition_curr = i + true_later;
        transitions = [transitions;transition_curr];
        
        % advance count past this
        i = transition_curr + n_periods;
        
    else
        
        % advance count  by 1
        i = i + 1;
        
    end
end

% remove transitions too close to the start or end of the record
transitions(transitions-time_to_take_spikes <=0) = [];
transitions(transitions+time_to_take_spikes > length(sleep)) = [];

% remove clusters if I want
ntransitions = length(transitions);
if rm_cluster
    trans_to_rm = zeros(ntransitions,1);
    all_trans_idx = (1:ntransitions)';
    for s = 1:ntransitions
        loo = all_trans_idx ~= s;
        if any(abs(transitions(s) - transitions(loo)) < time_to_take_spikes)
            trans_to_rm(s) = 1;
        end
    end
    transitions(logical(trans_to_rm)) = [];
    ntransitions = length(transitions);
end

% show it
days_per_index = 1/6/24;
times = linspace(0,length(sleep)*days_per_index,length(sleep));
if 0
    figure
    plot(times,sleep,'o')
    hold on
    for t = 1:length(transitions)
        plot([transitions(t)*days_per_index transitions(t)*days_per_index],ylim,'r--')
    end
end

bins = nan(length(transitions),time_to_take_spikes*2);
for t = 1:length(transitions)
    bins(t,:) = transitions(t) - time_to_take_spikes:transitions(t) + time_to_take_spikes-1;
end

end