function correlation_figure

%% Parameters
ss = 1; % just look at all sleep stages for simplicity
montage_names = {'car'};%{'machine','car','bipolar'};
all_montages = {'machine','car','bipolar'}; 

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
nmontages = length(montage_names);

% Establish network and univariate measures
networks = {'all_coh','all_pearson','all_plv','all_xcor','all_re'};
univariate = {'all_spikes','all_rl','all_bp','all_rel_bp','all_se', 'all_ll'};

%% Load data file
mt_data = load([inter_folder,'mt_out.mat']);
mt_data = mt_data.out;

%% Find non-empty pts
non_empty = find(cellfun(@(x) ~isempty(x), mt_data.all_bp(:,1,1)));
first_non_empty = non_empty(1);
npts = size(mt_data.all_bp,1);

%% Prep figure
figure

%% Subfigure A: show the correlation between different nodal measures (take average for connectivity measures)
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
            
            if nfreqi == 1
                net_namesi{icount} = sprintf('%s %s',strrep(net_names{in},'_',' '),montage_names{im});
            else
                net_namesi{icount} = sprintf('%s %s %s',strrep(net_names{in},'_',' '),montage_names{im},freq_names{f});
            end

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
                        if nfreqj == 1
                            net_namesj{jcount} = sprintf('%s %s',strrep(net_names{jn},'_',' '),montage_names{jm});
                        else
                            net_namesj{jcount} = sprintf('%s %s %s',strrep(net_names{jn},'_',' '),montage_names{jm},freq_names{jfreq});
                        end
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

% Remove those with all nans
% confirm first the only thing where this is true is relative bp broadband
% (which makes sense)
all_nans = isnan(nanmean(ic,2));
%assert(isequal(net_namesj(all_nans),{sprintf('rel bp %s broadband',montage_names{1})}))


nexttile
turn_nans_gray(ic(~all_nans,~all_nans))
xticks(1:size(ic,1)-1);
yticks(1:size(ic,1)-1);
%xticklabels(strrep(net_namesj(~all_nans),sprintf(' %s',montage_names{1}),''))
yticklabels(strrep(net_namesi(~all_nans),sprintf(' %s',montage_names{1}),''))
colorbar
clim([-1 1])
set(gca,'fontsize',15)
title(sprintf('Inter-feature correlation (all sleep stages, %s montage)',montage_names{1}))

%% Subfigure B: How much correlation is there across montages for different features?
thing_montages = nan(length(all_things),3);
thing_montages_sd = nan(length(all_things),3);
for in = 1:length(all_things)
    curr_net = mt_data.(all_things{in});
    temp_net = curr_net(:,1,1);
    if ismember(all_things{in},networks)
        nfreq = size(temp_net{first_non_empty},3); % how many frequencies
    else
        nfreq = size(temp_net{first_non_empty},2); % how many frequencies
    end

    net = cell(length(temp_net),3,3,nfreq); % npts,montage, ss

    % Loop over frequencies
    for f = 1:nfreq
    
        % Loop over montages
        for im = 1:3
            
            % Loop over sleep stages
            for is = 1:3
                cnet = curr_net(:,im,is);
                
                for ip = 1:length(cnet)
                
                    if isempty(cnet{ip}), continue; end
                    if ismember(all_things{in},networks)
                        net{ip,im,is,f} = (squeeze(nanmean(cnet{ip}(:,:,f),2))); % take the average connectivity for each electrode
                    else
                        net{ip,im,is,f} = (squeeze(cnet{ip}(:,f))); % just take the thing
                    end
    
                end
    
            end
    
        end

        


    end

    % Now do correlation across montages
    
    montage_corr = nan(3,npts,nfreq);
    which_montages = nan(3,2);
    for im = 1:3

        % will yield im=1->[2,3], im=2->[3,1], im=3->[1,2]
        m1 = mod(im,3)+1;
        m2 = mod(im+1,3)+1;
        which_montages(im,:) = [m1,m2];

        % Loop over frequencies and patients
        for f = 1:nfreq
            for ip = 1:npts
                mnet1 = net{ip,m1,ss,f};
                mnet2 = net{ip,m2,ss,f};
                if isempty(mnet1) || isempty(mnet2)
                    continue
                end
                montage_corr(im,ip,f) = corr(mnet1,mnet2,"rows","pairwise");
            end
        end

        
    end

    % Average the correlations across frequencies and patients
    for im = 1:3
        thing_montages(in,im) = squeeze(nanmean(montage_corr(im,:,:),[2 3]));
        thing_montages_sd(in,im) = squeeze(nanmean(nanstd(montage_corr(im,:,:),[],2),3));
    end


end

nexttile
cols = colormap(gca,lines(3));
lp = nan(3,1);
offset = 0.1;
for i = 1:size(thing_montages,1)
    for j = 1:size(thing_montages,2)
        lp(j) = errorbar(i+offset*j-offset*2,thing_montages(i,j),thing_montages_sd(i,j),...
            'o','color',cols(j,:),'linewidth',2,'markersize',12);
        hold on

    end
end
xticks(1:size(thing_montages,1))
xticklabels(net_names)
ylabel('Correlation (r)')
labels = cell(3,1);
for i =1:3
    labels{i} = sprintf('%s-%s',strrep(all_montages{which_montages(i,1)},'_',' '),strrep(all_montages{which_montages(i,2)},'_',' '));
end
legend(lp,labels,'location','southeast','fontsize',15)
set(gca,'fontsize',15)
title('Inter-reference feature correlation')


%% Subfigure C: How much correlaiton is there across sleep stages for different features?


end