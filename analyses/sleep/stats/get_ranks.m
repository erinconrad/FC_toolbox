function [all_ranks,all_soz_ranks,nchance,all,successes] = ...
    get_ranks(things,sozs,rates,which,min_rate)

%% Initialize data things
npts = length(things);
all_ranks = cell(npts,1);
all_soz_ranks = cell(npts,1);
nchance = nan(npts,1);
all = nan(npts,1);
successes = nan(npts,1);

for i = 1:npts
    curr_things = things{i};
    curr_soz = sozs{i};
    curr_rates = rates{i};
    
    % remove nans
    nan_things = isnan(curr_things) | isnan(curr_rates);
    curr_things(nan_things) = [];
    curr_soz(nan_things) = [];
    curr_soz = logical(curr_soz);
    
    curr_rates(nan_things) = [];
    
    % get rank of soz
    if strcmp(which,'rate')

        % sort the rates in descending order
        [curr_things,I] = sort(curr_things,'descend');
    elseif strcmp(which,'rl')

        % also remove those with few spikes. I do this because I think
        % the RL for electrodes with few spikes is unreliable
        spikey = curr_rates > min_rate;
        curr_things = curr_things(spikey);
        curr_soz = curr_soz(spikey);

        % sort the RL in ascending order
        [curr_things,I] = sort(curr_things,'ascend');
    end
    ranks = 1:length(curr_things);
    ranks(I) = ranks;
    soz_ranks = ranks(curr_soz);
    not_soz_ranks = ranks(~curr_soz);
    nchance(i) = length(ranks)/2;
    all(i) = nanmedian(soz_ranks);
    
    successes(i) = ranking_binomial_test(ranks,curr_soz);
    
    %{
    plot(i+0*randn(length(not_soz_ranks),1),not_soz_ranks,'o',...
        'color',grayColor,'markersize',markersize);
    hold on
    plot(i+0*randn(length(soz_ranks),1),soz_ranks,'o',...
        'color',myColours(1,:),'linewidth',2,'markersize',markersize);
    %}
    all_ranks{i} = ranks;
    all_soz_ranks{i} = soz_ranks;
    
end

end