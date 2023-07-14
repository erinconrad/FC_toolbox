function correlation_figure

%% Parameters
ss = 3; % just look at sleep
montage_names = {'car'};%{'machine','car','bipolar'};
all_montages = {'machine','car','bipolar'}; 
which_pts = 'hup';

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
univariate = {'all_spikes','all_rl','all_bp','all_se', 'all_ll'};

%% Initialize results file
fname = [plot_folder,'supplemental_results.html'];
fid = fopen(fname,'a');
fprintf(fid,'<p><br><u><i>Many interictal EEG features are highly correlated, and the choice of reference affects features</i></u></br>');
fprintf(fid,['We examined the correlation between interictal EEG features. ']);

%% Load data file
mt_data = load([inter_folder,'mt_out_epilepsy_laterality.mat']);
mt_data = mt_data.out;

%% Get appropriate patients
switch which_pts
    case 'hup'
        pts = contains(mt_data.all_names,'HUP');
    case 'musc'
        pts = contains(mt_data.all_names,'MP');
end
npts = length(pts);

%% Find non-empty pts
non_empty = find(cellfun(@(x) ~isempty(x), mt_data.all_bp(pts,1,1)));
first_non_empty = non_empty(1);
%npts = size(mt_data.all_bp,1);

%% Prep figure
figure
set(gcf,'position',[10 10 1300 1200])
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')

%% Subfigure A: show the correlation between different nodal measures (take average for connectivity measures)
all_things = [networks,univariate];
net_names = cellfun(@(x) strrep(x,'all_',''),all_things,'UniformOutput',false);

% get number of networks
nnet = 0;
for in = 1:length(all_things)
    curr_net = mt_data.(all_things{in});
    curr_net = curr_net(pts,1,ss);
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
        curr_neti = curr_neti(pts,im,ss);

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
                    curr_netj = curr_netj(pts,jm,ss);

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

clean_net_names = strrep(net_namesi(~all_nans),sprintf(' %s',montage_names{1}),'');
nexttile
turn_nans_gray(ic(~all_nans,~all_nans))
xticks(1:size(ic,1));
yticks(1:size(ic,1));
%xticklabels(strrep(net_namesj(~all_nans),sprintf(' %s',montage_names{1}),''))
xticklabels(cellfun(@greek_letters_plots,clean_net_names,'uniformoutput',false))
yticklabels(cellfun(@greek_letters_plots,clean_net_names,'uniformoutput',false))
colorbar
clim([-1 1])
set(gca,'fontsize',15)
title(sprintf('Inter-feature correlation (electrode contact level)'))

fprintf(fid,['We first measured inter-feature correlation on an electrode contact-level. To do '...
    'this, we converted bivariate features to electrode contact-specific '...
    'univariate features by taking the average edge weight across all other electrode contacts. (For '...
    'this analysis, we did not restrict the average to contacts only on the same electrode). This yielded '...
    'a single measure for each patient, feature, electrode contact, and reference. '...
    'We then calculated the Pearson correlation across electrode contacts between all features and choices of reference, '...
    'yielding a <i>N</i><sub>features x references</sub> x <i>N</i><sub>features x references</sub> '...
    'correlation matrix for each patient, where <i>N</i><sub>features x references</sub> is the '...
    'number of features times the number of references (30 x 3 = 90). Fig. S1A shows the '...
    'average inter-feature correlation matrix across patients for a single choice of reference (common average). '...
    'There were often high correlations between different frequency band measurements of the same feature. '... ...
    'There were also often high correlations or anti-correlations between different features, such as between '...
    'coherence and phase-locking value, and between coherence and relative entropy (negative correlation).'...
    ' Spike rates were moderately positively correlated with relative entropy and bandpower.']);

%% Subfigure B: How much correlation is there across montages for different features?
which_freq = 6; % broadband
thing_montages = nan(length(all_things),3);
thing_montages_sd = nan(length(all_things),3);
for in = 1:length(all_things)
    curr_net = mt_data.(all_things{in});
    temp_net = curr_net(pts,1,1);
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

    % Average the correlations across patients, pick broadband frequency
    for im = 1:3
        if nfreq == 1
            thing_montages(in,im) = squeeze(nanmean(montage_corr(im,:,1),2));
            thing_montages_sd(in,im) = squeeze((nanstd(montage_corr(im,:,1),[],2)));
        else
            thing_montages(in,im) = squeeze(nanmean(montage_corr(im,:,which_freq),2));
            thing_montages_sd(in,im) = squeeze((nanstd(montage_corr(im,:,which_freq),[],2)));
        end
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

ylim([-1 1.3])
plot(xlim,[0 0],'k--')
xticks(1:size(thing_montages,1))
xlim([0.5 size(thing_montages,1)+0.5])
xticklabels(cellfun(@greek_letters_plots, cellfun(@(x) strrep(x,'_',' '),net_names,'uniformoutput',false),'uniformoutput',false))
ylabel('Correlation (r)')
labels = cell(3,1);
for i =1:3
    labels{i} = sprintf('%s-%s',strrep(all_montages{which_montages(i,1)},'_',' '),strrep(all_montages{which_montages(i,2)},'_',' '));
end
set(gca,'fontsize',15)
for i = 1:size(thing_montages,1)
    if i<size(thing_montages,1)
        plot([i+0.5 i+0.5],ylim,'k--')
    end
end

legend(lp,labels,'location','southeast','fontsize',15)

title('Inter-reference feature correlation (contact level)')

fprintf(fid,[' Fig. S1B shows the mean (standard deviation) inter-reference feature correlation '...
    'across patients for different features. Only the broadband frequency-measured features '...
    'were considered for this analysis. Correlations varied across features, and the correlation '...
    'in features between machine reference and common average reference tended to be higher than that between '...
    'bipolar reference and either other reference. Several correlations were less than 0.5, '...
    'suggesting that the choice of reference strongly influences interictal EEG feature calculations, and may '...
    'have a stronger influence than the choice of connectivity measurement.</p>']);

%% Now do lr_mt to get AI features
which_sleep_stages = 3;
[T,features] =  lr_mt(which_sleep_stages);

% Restrict to correct hospital
switch which_pts
    case 'hup'
        T = T(contains(T.names,'HUP'),:);
    case 'musc'
        T = T(contains(T.names,'MP'),:);
end

%% Subfigure C: Inter-AI correlation
% Restrict to desired montage and no SD
montage_features = contains(features,sprintf('%s ',montage_names{1}));
fT = T(:,features);
fT = fT(:,montage_features);
no_sd = ~contains(features(montage_features),'SD');
fT = fT(:,no_sd);
no_sd_names = features(montage_features); no_sd_names = no_sd_names(no_sd);

% clean the names
no_sd_names = strrep(no_sd_names,sprintf(' %s',montage_names{1}),'');
no_sd_names = strrep(no_sd_names,' sleep','');

% Do inter-feature correlation
aic = corr(table2array(fT),table2array(fT),'rows','pairwise');

% re-order things so the names match with the first subfigure
[lia,locb] = ismember(clean_net_names,no_sd_names);
assert(all(lia==1))
no_sd_names = no_sd_names(locb);
aic = aic(locb,locb);

% plot it
nexttile
turn_nans_gray(aic)
xticks(1:size(aic,1));
yticks(1:size(aic,1));
%xticklabels(strrep(net_namesj(~all_nans),sprintf(' %s',montage_names{1}),''))
xticklabels(cellfun(@greek_letters_plots,no_sd_names,'uniformoutput',false))
yticklabels(cellfun(@greek_letters_plots,no_sd_names,'uniformoutput',false))
colorbar
clim([-1 1])
set(gca,'fontsize',15)
title(sprintf('Inter-feature asymmetry index correlation (patient level)'))



%% Subfigure D: Inter-reference AI correlation
if which_freq == 6
    freq_name = 'broadband';
else
    error('why are you messing with frequency')
end
ai_montages = nan(length(all_things),3);
which_montages = nan(3,2);
all_feature_names = strrep(all_things,'all_','');
for in = 1:length(all_things)

    curr_features = features(contains(features,sprintf('%s ',all_feature_names{in})));

    % remove SD
    curr_features(contains(curr_features,'SD')) = [];

    % Get only desired frequency -  check if there are more than 3, that
    % means there are multiple frequencies
    if length(curr_features) == 18
        curr_features(~contains(curr_features,freq_name)) = [];
    elseif length(curr_features) == 3 % do nothign
    else
        error('what')
    end
    

    for im = 1:3
    
        % will yield im=1->[2,3], im=2->[3,1], im=3->[1,2]
        m1 = mod(im,3)+1;
        m2 = mod(im+1,3)+1;
        which_montages(im,:) = [m1,m2];
    
        mname1 = all_montages{m1};
        mname2 = all_montages{m2};

        mname1_feature = curr_features(contains(curr_features,mname1));
        mname2_feature = curr_features(contains(curr_features,mname2));

        % correlate the AIs
        ai_montages(in,im) = corr(T.(mname1_feature{1}),T.(mname2_feature{1}),'rows','pairwise');

    
    end
end

% Ensure names align
assert(isequal(all_feature_names,cellfun(@(x) strrep(x,'_',' '),net_names,'uniformoutput',false)))

% Plot it
nexttile
cols = colormap(gca,lines(3));
lp = nan(3,1);
offset = 0.1;
for i = 1:size(ai_montages,1)
    
    for j = 1:size(ai_montages,2)
        lp(j) = plot(i+offset*j-offset*2,ai_montages(i,j),...
            'o','color',cols(j,:),'linewidth',2,'markersize',12);
        hold on

    end

    
end
ylim([-1 1])
plot(xlim,[0 0],'k--')
set(gca,'fontsize',15)
xticks(1:size(ai_montages,1))
xlim([0.5 size(ai_montages,1)+0.5])
xticklabels(all_feature_names)
ylabel('Correlation (r)')
labels = cell(3,1);
for i =1:3
    labels{i} = sprintf('%s-%s',strrep(all_montages{which_montages(i,1)},'_',' '),strrep(all_montages{which_montages(i,2)},'_',' '));
end
for i = 1:size(ai_montages,1)
    if i<size(ai_montages,1)
        plot([i+0.5 i+0.5],ylim,'k--')
    end
end


legend(lp,labels,'location','southeast','fontsize',15)

title('Inter-reference asymmetry index correlation (patient level)')

fprintf(fid,['<p>We repeated this inter-feature correlation analysis, this time '...
    'measuring the correlation between the AI of the different features. '...
    'In this analysis, given that there was a single AI measurement for each '...
    '<i>patient</i> (rather than each electrode contact), inter-feature Pearson correlations '...
    'were calculated across patients, rather than across electrode contacts. '...
    'Similar trends were observed as in the electrode contact-level analysis. '...
    'Fig. S1C shows the inter-feature AI correlation, again restricting analysis to common '...
    'average reference. Again, we observed high inter-feature correlation for several features. '...
    'Fig. S1D shows the inter-reference AI correlation (there are no error bars because '...
    ' there is only a single correlation value across all patients for this analysis). '...
    'Again, we observed high variability of inter-reference correlations '...
    'across features.</p>']);

%% Add subtitles
annotation('textbox',[0 0.905 0.1 0.1],'String','A','LineStyle','none','fontsize',20)
annotation('textbox',[0.53 0.905 0.1 0.1],'String','B','LineStyle','none','fontsize',20)
annotation('textbox',[0 0.41 0.1 0.1],'String','C','LineStyle','none','fontsize',20)
annotation('textbox',[0.53 0.41 0.1 0.1],'String','D','LineStyle','none','fontsize',20)

print(gcf,[plot_folder,'FigS1'],'-dpng')


end