function [transitions,bins] = designate_histogram(sleep,wake,n_periods,min_same,...
    later_search,time_to_take_spikes,rm_cluster,ad_norm,disc,name,out_folder)

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980];

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
    if sum(wake(before) == 1)/n_periods > min_same && ...
            sum(sleep(after) == 1)/n_periods > min_same
        
        sum_appropriate = zeros(later_search+1,1);
        % this is a potential transition point, search later to see if
        % there is a better one
        for later = 0:later_search % up to 120 minutes later
            
            % count matching sleep/wake for this transition point
            new_before = i+later - n_periods:i+later-1;
            new_after = i+later:i+later+n_periods;
            
            sum_appropriate(later+1) = sum(wake(new_before) == 1) + ...
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
if 1
    figure
    set(gcf,'position',[10 10 1400 600])
    tt = tiledlayout(2,1,'padding','tight','tilespacing','tight');
    
    nexttile
    plot(times,ad_norm,'linewidth',2)
    hold on
    plot(xlim,[disc disc],'k--','linewidth',2)
    ylabel('Normalized alpha delta ratio')
    xlabel('Days')
    title('Alpha delta ratio')
    
    set(gca,'fontsize',15)
    
    nexttile
    plot(times(wake==1),0,'o','markersize',4,'color',colors(1,:))
    hold on
    
    plot(times(sleep==1),1,'o','markersize',4,'color',colors(2,:))
    yticks([0 1])
    ylim([-1 2])
    yticklabels({'Wake','Sleep'})
    hold on
    for t = 1:length(transitions)
        plot([transitions(t)*days_per_index transitions(t)*days_per_index],ylim,...
            'k--','linewidth',2)
    end
    set(gca,'fontsize',15)
    xlabel('Days')
    title('Wake-sleep transition times')
    
    title(tt,name,'fontsize',20,'fontweight','bold')
    
    print([out_folder,name],'-dpng');
    close(gcf)
end

bins = nan(length(transitions),time_to_take_spikes*2);
for t = 1:length(transitions)
    bins(t,:) = transitions(t) - time_to_take_spikes:transitions(t) + time_to_take_spikes-1;
end

end