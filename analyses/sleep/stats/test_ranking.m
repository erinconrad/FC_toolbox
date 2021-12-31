function mcout = test_ranking(things,soz,nb,which,rates,min_rate)

% Need to test that I am indexing correctly by grabbing electrode names

npts = length(soz);


%% Get MC median rankings (bootstrap and true)
median_ranking_all = nan(nb+1,1); % first one is true
all_rankings = nan(npts,1);
for ib = 0:nb
    if mod(ib,100) == 0
        fprintf('\nDoing %d out of %d\n',ib,nb);
    end
    temp_all_rankings = nan(npts,1);
    no_soz = zeros(npts,1);
    % Loop over patients
    for ip = 1:npts
    
        % Sozs for that patient
        curr_soz = soz{ip};
        chs = 1:length(things{ip});
        curr_soz_bin = curr_soz;
        curr_soz_bin(ismember(chs,curr_soz)) = 1;
        
        % remove from thing and soz cases of nan values
        curr_things = things{ip};
        nan_things = isnan(curr_things);
        curr_things(nan_things) = [];
        curr_rates = rates{ip};
        curr_rates(nan_things) = [];
        
        curr_soz_bin(nan_things) = [];
        curr_soz_bin = logical(curr_soz_bin);
        
        if strcmp(which,'rate')

            % sort the rates in descending order
            [curr_things,I] = sort(curr_things,'descend');
        elseif strcmp(which,'rl')
            
            % also remove those with few spikes. I do this because I think
            % the RL for electrodes with few spikes is unreliable
            spikey = curr_rates > min_rate;
            curr_things = curr_things(spikey);
            curr_soz_bin = curr_soz_bin(spikey);
            
            % sort the RL in ascending order
            [curr_things,I] = sort(curr_things,'ascend');
        end

        % Get number of remaining soz electrodes and remaining all
        % electrodes
        nsoz = sum(curr_soz_bin);
        nelecs = length(curr_things);
        if nsoz == 0
            no_soz(ip) = 1;
        end
    
        % Turn I into ranks
        ranks = 1:length(curr_things);
        ranks(I) = ranks;

        % if ib == 0, then it's the true thing. If not, it's a Monte Carlo
        % iteration
        if ib ~= 0
            % choose a random sample (without replacement) of electrodes equal
            % in number to nsoz
            fake_soz = randsample(nelecs,nsoz);
        else
            fake_soz = find(curr_soz_bin);
            
        end
        
        % get ranks of SOZ electrodes
        soz_ranks = median(ranks(fake_soz));
        temp_all_rankings(ip) = soz_ranks;
        
        if ib == 0
            all_rankings(ip) = soz_ranks;
        end
    end

    % remove patients who have no SOZ
    no_soz = logical(no_soz);
    median_ranking_all(ib+1) = median(temp_all_rankings(~no_soz));
    
end

median_ranking_true = median_ranking_all(1);
median_ranking_mc = median_ranking_all(2:end);
pval = (sum(median_ranking_mc <= median_ranking_true)+1)/(nb+1);

if 0
    plot(sort(median_ranking_mc),'o')
    hold on
    plot(xlim,[median_ranking_true median_ranking_true],'linewidth',2)
    title(get_p_text(pval))
    legend({'Monte Carlo','True'})
end

mcout.pval = pval;
mcout.all_rankings = all_rankings;
mcout.median_ranking_mc = median_ranking_mc;
mcout.median_ranking_true = median_ranking_true;


end