function [T,features,way,dur,sample,ss,durations] =  lr_mt_multitime(mt_data,which_sleep_stages)

%% Parameters
only_hup = 0;
which_montages = [1 2 3];
do_little_plots = 0;


all_missing = cellfun(@isempty,mt_data.all_spikes(:,1,1));
names = mt_data.all_names;
npts = length(names);
surgery = mt_data.all_surgery;
resection_lat = mt_data.all_resec_lat;
resection_loc = mt_data.all_resec_loc;
ablation_lat = mt_data.all_ablate_lat;
ablation_loc = mt_data.all_ablate_loc;
engel_yr1 = mt_data.all_engel(:,1);
engel_yr2 = mt_data.all_engel(:,2);
ilae_yr1 = mt_data.all_ilae(:,1);
ilae_yr2 = mt_data.all_ilae(:,2);
soz_lats = mt_data.all_soz_lat;
soz_locs = mt_data.all_soz_loc;
disconnected = mt_data.all_disconnected;
all_n_wake_sleep_connected = mt_data.all_n_wake_sleep_connected;

all_spikes = mt_data.spikes_subsample;
durations = mt_data.durations;

% Find and exclude patients for whom bulk of record is disconnected
most_disconnected = sum(disconnected == 1,2) >= 0.9* size(disconnected,2);

% Identify patients for whom there is little wake or sleep connected
no_wake = all_n_wake_sleep_connected(:,1) < 5;
no_sleep = all_n_wake_sleep_connected(:,2) < 5;
n_wake = all_n_wake_sleep_connected(:,1);
n_sleep = all_n_wake_sleep_connected(:,2);
n_connected = sum(disconnected == 0,2);

% get number of symmetric labels
ref_labels = mt_data.all_labels(:,1);
n_symmetric = cellfun(@length,ref_labels);


%% Clean SOZ localizations and lateralities
soz_lats(cellfun(@isempty,soz_lats)) = {''};
soz_locs(cellfun(@isempty,soz_locs)) = {''};
soz_lats(strcmp(soz_lats,'diffuse')) = {'bilateral'}; % make diffuse be the same as bilateral
soz_locs(contains(soz_locs,'temporal')) = {'temporal'}; % make any temporal be temporal

%% Consensus ablation or resection lat
surg_lat = cell(npts,1);
for i = 1:npts
    if isempty(resection_lat{i}) && isempty(ablation_lat{i})
    elseif strcmp(resection_lat{i},'na') && strcmp(ablation_lat{i},'na')
    elseif strcmp(resection_lat{i},'na')
        surg_lat{i} = ablation_lat{i};
    elseif strcmp(ablation_lat{i},'na')
        surg_lat{i} = resection_lat{i};
    else
        if ~strcmp(resection_lat{i},ablation_lat{i})
            error('what');
        end
        surg_lat{i} = resection_lat{i};
    end

end

%% Consensus ablation or reseciton loc
surg_loc = cell(npts,1);
for i = 1:npts
    if isempty(resection_loc{i}) && isempty(ablation_loc{i})
    elseif strcmp(resection_loc{i},'na') && strcmp(ablation_loc{i},'NA')
    elseif strcmp(resection_loc{i},'ATL')
        surg_loc{i} = 'temporal';
    elseif contains(ablation_loc{i},'temporal')
        surg_loc{i} = 'temporal';
    else
        surg_loc{i} = 'other';
    end
end


%% Get features
% Initialize table
Ts = table(names,engel_yr1,engel_yr2,ilae_yr1,ilae_yr2,surgery,surg_lat,surg_loc,soz_locs,soz_lats,no_wake,no_sleep,n_wake,n_sleep,n_connected,n_symmetric);
features = {};
way = [];
dur = [];
sample = [];
ss = [];

for which_sleep_stage = which_sleep_stages% all = 1, wake =2, sleep = 3;
    sleep_text = sprintf('_ss%d',which_sleep_stage);
    

    for which_montage =which_montages % machine = 1,car = 2, bipolar = 3
        
        if which_montage == 1
            montage_text = 'machine';
        elseif which_montage == 2
            montage_text = 'car';
        elseif which_montage == 3
            montage_text = 'bipolar';
        end
        labels = mt_data.all_labels(:,which_montage);
        
       
        spikes = all_spikes(:,which_montage,which_sleep_stage,:,:,:);
        ndurations = size(spikes,5);
        nsamples = size(spikes,6);

        % Loop over the two ways of doing subsampling
        for iw = 1:2

            way_text = sprintf('_way%d',iw);

            % Loop over durations
            for id = 1:ndurations

                dur_text = sprintf('_dur%d',id);
                
                % Loop over subsamplings
                for is = 1:nsamples

                    samp_text = sprintf('_samp%d_',is);
                    
                    % Get specific spikes
                    thing = squeeze(spikes(:,1,1,iw,id,is));
                    which_thing = {'spikes'};
                    uni = 1;
                    last_dim = 1;

                    % Get AI
                    ai = cell2mat(cellfun(@(x,y,z,w,a,d) ...
                    calc_ai_ns(x,y,z,w,a,d,uni,last_dim,which_thing,subplot_path,do_little_plots),...
                    labels,thing,names,mt_data.all_labels(:,1),atropos,dkt,'uniformoutput',false));

                    % Fill up table
                    feat_name = [which_thing{1},'_',montage_text,'_',sleep_text,...
                        way_text,dur_text,samp_text];
                    features = [features,feat_name];
                    way = [way;iw];
                    dur = [dur;id];
                    sample = [sample;is];
                    ss = [ss;which_sleep_stage];

                    Ts = addvars(Ts,ai);
                    Ts = splitvars(Ts,'ai','NewVariableNames',feat_name);
                end

            end

        end
              
  
    end

end


%% Prep  output table
% First, remove those rows missing all columns
T = Ts(~all_missing & ~most_disconnected,:);



end