function network_correlation

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
subplot_path = [plot_folder,'ai_subplots/'];
if ~exist(subplot_path,'dir')
    mkdir(subplot_path)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% Frequency band names
freq_names = {'delta','theta','alpha','beta','gamma','broadband'};
montage_names = {'machine','car','bipolar'}; nmontages = length(montage_names);

%% Load data file
mt_data = load([inter_folder,'mt_out.mat']);
mt_data = mt_data.out;

%% Find non-empty pts
non_empty = find(cellfun(@(x) ~isempty(x), mt_data.all_bp(:,1,1)));
first_non_empty = non_empty(1);


%% Initialize figure
%figure

%% Network correlation
ss = 1; % all sleep stages
networks = {'all_coh','all_pearson','all_plv','all_xcor','all_re'};
net_names = cellfun(@(x) strrep(x,'all_',''),networks,'UniformOutput',false);

% get number of networks
nnet = 0;
for in = 1:length(networks)
    curr_net = mt_data.(networks{in});
    curr_net = curr_net(:,1,ss);
    nnet = nnet + nmontages*size(curr_net{first_non_empty},3);
end

% prep inter-network correlation matrix
ic = nan(nnet,nnet);

icount = 0;
net_namesi = cell(nnet,1);
net_namesj = cell(nnet,1);
% loop over montages
for im = 1:nmontages

    % Loop over networks
    for in = 1:length(networks)

        curr_neti = mt_data.(networks{in});
        curr_neti = curr_neti(:,im,ss);
        nfreqi = size(curr_neti{first_non_empty},3); % how many frequencies

        for f = 1:nfreqi
            icount = icount + 1;


            neti = cell(length(curr_neti),1);
            for ip = 1:length(curr_neti)
                if isempty(curr_neti{ip}), continue; end
                neti{ip} = (squeeze(curr_neti{ip}(:,:,f)));
            end
            
            net_namesi{icount} = sprintf('%s %s %s',net_names{in},montage_names{im},freq_names{f});

            jcount = 0;
            for jm = 1:nmontages
                for jn = 1:length(networks)
                    curr_netj = mt_data.(networks{jn});
                    curr_netj = curr_netj(:,jm,ss);
                    nfreqj = size(curr_netj{first_non_empty},3); % how many frequencies

                    for jfreq = 1:nfreqj
                        jcount =jcount +1;
                        netj= cell(length(curr_neti),1);
                        for ip = 1:length(curr_netj)
                            if isempty(curr_netj{ip}), continue; end
                            netj{ip} = (squeeze(curr_netj{ip}(:,:,jfreq)));
                        end
                        net_namesj{jcount} = sprintf('%s %s %s',net_names{jn},montage_names{jm},freq_names{jfreq});
                        all_pts_corr = nan(length(curr_netj),1);
                        for ip = 1:length(curr_netj)
                            
                            if isempty(netj{ip}), continue; end
                            all_pts_corr(ip) = corr(wrap_or_unwrap_adjacency_fc_toolbox(neti{ip}),...
                                wrap_or_unwrap_adjacency_fc_toolbox(netj{ip}),'rows','pairwise');
                        end
                        ic(icount,jcount) = nanmean(all_pts_corr);
                    end

                end
            end

        end


    end

end

figure
turn_nans_gray(ic)
xticks(1:size(ic,1));
yticks(1:size(ic,1));
xticklabels(net_namesj)
yticklabels(net_namesi)
colorbar
clim([-1 1])
set(gca,'fontsize',15)


%% Univariate measures
univariate = {'all_spikes','all_rl','all_bp','all_rel_bp','all_se', 'all_ll'};
net_names = cellfun(@(x) strrep(x,'all_',''),univariate,'UniformOutput',false);

% get number of networks
nnet = 0;
for in = 1:length(univariate)
    curr_net = mt_data.(univariate{in});
    curr_net = curr_net(:,1,ss);
    nnet = nnet + nmontages*size(curr_net{first_non_empty},2);
end

% prep inter-network correlation matrix
ic = nan(nnet,nnet);

icount = 0;
net_namesi = cell(nnet,1);
net_namesj = cell(nnet,1);
% loop over montages
for im = 1:nmontages

    % Loop over networks
    for in = 1:length(univariate)

        curr_neti = mt_data.(univariate{in});
        curr_neti = curr_neti(:,im,ss);
        nfreqi = size(curr_neti{first_non_empty},2); % how many frequencies

        for f = 1:nfreqi
            icount = icount + 1;


            neti = cell(length(curr_neti),1);
            for ip = 1:length(curr_neti)
                if isempty(curr_neti{ip}), continue; end
                neti{ip} = (squeeze(curr_neti{ip}(:,f)));
            end
            
            net_namesi{icount} = sprintf('%s %s %s',net_names{in},montage_names{im},freq_names{f});

            jcount = 0;
            for jm = 1:nmontages
                for jn = 1:length(univariate)
                    curr_netj = mt_data.(univariate{jn});
                    curr_netj = curr_netj(:,jm,ss);
                    nfreqj = size(curr_netj{first_non_empty},2); % how many frequencies

                    for jfreq = 1:nfreqj
                        jcount =jcount +1;
                        netj= cell(length(curr_neti),1);
                        for ip = 1:length(curr_netj)
                            if isempty(curr_netj{ip}), continue; end
                            netj{ip} = (squeeze(curr_netj{ip}(:,jfreq)));
                        end
                        net_namesj{jcount} = sprintf('%s %s %s',net_names{jn},montage_names{jm},freq_names{jfreq});
                        all_pts_corr = nan(length(curr_netj),1);
                        for ip = 1:length(curr_netj)
                            
                            if isempty(netj{ip}), continue; end
                            all_pts_corr(ip) = corr((neti{ip}),...
                                (netj{ip}),'rows','pairwise');
                        end
                        ic(icount,jcount) = nanmean(all_pts_corr);
                    end

                end
            end

        end


    end

end

figure
turn_nans_gray(ic)
xticks(1:size(ic,1));
yticks(1:size(ic,1));
xticklabels(net_namesj)
yticklabels(net_namesi)
colorbar
clim([-1 1])
set(gca,'fontsize',15)

%% Nodal measures (Everything!)
all_things = [networks,univariate];
net_names = cellfun(@(x) strrep(x,'all_',''),all_things,'UniformOutput',false);

% get number of networks
nnet = 0;
for in = 1:length(all_things)
    curr_net = mt_data.(all_things{in});
    curr_net = curr_net(:,1,ss);
    if ismember(all_things{in},networks)
        nnet = nnet + nmontages*size(curr_net{first_non_empty},3);
    else
        nnet = nnet + nmontages*size(curr_net{first_non_empty},2);
    end
end

% prep inter-network correlation matrix
ic = nan(nnet,nnet);

icount = 0;
net_namesi = cell(nnet,1);
net_namesj = cell(nnet,1);
% loop over montages
for im = 1:nmontages

    % Loop over networks
    for in = 1:length(all_things)

        curr_neti = mt_data.(all_things{in});
        curr_neti = curr_neti(:,im,ss);

        if ismember(all_things{in},networks)
            nfreqi = size(curr_neti{first_non_empty},3); % how many frequencies
        else
            nfreqi = size(curr_neti{first_non_empty},2); % how many frequencies
        end

        for f = 1:nfreqi
            icount = icount + 1;


            neti = cell(length(curr_neti),1);
            for ip = 1:length(curr_neti)
                if isempty(curr_neti{ip}), continue; end
                if ismember(all_things{in},networks)
                    neti{ip} = (squeeze(nanmean(curr_neti{ip}(:,:,f),2))); % take the average connectivity for each electrode
                else
                    neti{ip} = (squeeze(curr_neti{ip}(:,f))); % just take the thing
                end
                    
            end
            
            net_namesi{icount} = sprintf('%s %s %s',net_names{in},montage_names{im},freq_names{f});

            jcount = 0;
            for jm = 1:nmontages
                for jn = 1:length(all_things)
                    curr_netj = mt_data.(all_things{jn});
                    curr_netj = curr_netj(:,jm,ss);

                    if ismember(all_things{jn},networks)
                        nfreqj = size(curr_netj{first_non_empty},3); % how many frequencies
                    else
                        nfreqj = size(curr_netj{first_non_empty},2); % how many frequencies
                    end

                    for jfreq = 1:nfreqj
                        jcount =jcount +1;
                        netj= cell(length(curr_neti),1);
                        for ip = 1:length(curr_netj)
                            if isempty(curr_netj{ip}), continue; end
                            if ismember(all_things{jn},networks)
                                netj{ip} = (squeeze(nanmean(curr_netj{ip}(:,:,jfreq),2)));
                            else
                                netj{ip} = (squeeze(curr_netj{ip}(:,jfreq)));
                            end
                        end
                        net_namesj{jcount} = sprintf('%s %s %s',net_names{jn},montage_names{jm},freq_names{jfreq});
                        all_pts_corr = nan(length(curr_netj),1);
                        for ip = 1:length(curr_netj)
                            
                            if isempty(netj{ip}), continue; end
                            all_pts_corr(ip) = corr((neti{ip}),...
                                (netj{ip}),'rows','pairwise');
                        end
                        ic(icount,jcount) = nanmean(all_pts_corr);
                    end

                end
            end

        end


    end

end

figure
turn_nans_gray(ic)
xticks(1:size(ic,1));
yticks(1:size(ic,1));
xticklabels(net_namesj)
yticklabels(net_namesi)
colorbar
clim([-1 1])
set(gca,'fontsize',15)
end