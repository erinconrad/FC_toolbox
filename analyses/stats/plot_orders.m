function plot_orders(things,sozs,rates,which,min_rate)

myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];
grayColor = 0.8*[1 1 1];
markersize = 5;

npts = length(things);
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
    
    plot(i+0*randn(length(not_soz_ranks),1),not_soz_ranks,'o',...
        'color',grayColor,'markersize',markersize);
    hold on
    plot(i+0*randn(length(soz_ranks),1),soz_ranks,'o',...
        'color',myColours(1,:),'linewidth',2,'markersize',markersize);
    
    
end


end