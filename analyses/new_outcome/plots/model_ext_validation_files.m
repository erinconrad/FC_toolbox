function model_ext_validation_files

%% Parameters
pca_perc = 95; % the percent variance to explain for pca
rm_non_temporal = 1; % remove patients who are not temporal
rm_wake = 1; % don't include wake segments
which_refs = {'machine','car','bipolar'};

locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

for ir = 1:length(which_refs)

    file_name = sprintf('ext_models_%s.mat',which_refs{ir});

    %% Main analysis and good spike analysis
    for ia = 1:2
    
        if ia == 1 % all spikes
            approach(ia).type = 'all spikes';
            % Run the lr_mt to extract AI features
            if rm_wake == 1
                [T,features] =  lr_mt(3); % the 3 refers to only looking at sleep
            else
                error('why are you doing this?')
            end
    
        else % good spikes only
            approach(ia).type = 'good spikes, PPV above 0.5';
            % Run the lr_mt to extract AI features, but restrict to good spikes
            % only
            if rm_wake == 1
                [T,features] =  lr_mt(3,1); % the 3 refers to only looking at sleep
            else
                error('why are you doing this?')
            end
    
        end
    
        
        % Remove those without a response (soz_lats is the response variable)
        empty_class = cellfun(@isempty,T.soz_lats);
        T(empty_class,:) = [];
        
        % Establish HUP and MUSC as training and testing, respectively
        train  = contains(T.names,'HUP');
        test  = contains(T.names,'MP');
        
        % Establish model types
        model(1).type = 'All features';
        model(2).type = 'Spikes';
        model(3).type = 'Binary spikes';
        nmodels = length(model);
        
        for im = 1:nmodels
            model(im).val(1).description = 'HUP cross-validation';
            model(im).val(2).description = 'MUSC external validation';
        end
        
        for im = 1:nmodels
            for iv = 1:2
                model(im).val(iv).side(1).description = 'Left vs right/bilateral';
                model(im).val(iv).side(2).description = 'Right vs left/bilateral';
            end
        end
        
        % Run all the models
        fprintf('\nDoing main models...');
        tic
        % Loop over model types
        for im = 1:nmodels
            just_spikes = im - 1; % if full model, just_spikes = 0, if spikes, just_spikes = 1, if binary, just_spikes = 2
            
           
            % Loop over internal vs external validation
            for iv = 1:2
        
                % Loop over sides
                for is = 1:2
        
                    % run the internal cross validation
                    if iv == 1
                        out = classifier_wrapper(T(train,:),features,pca_perc,is,...
                            just_spikes,rm_non_temporal,[],which_refs{ir});
            
                    % run the external test
                    elseif iv == 2
                        out = validation_classifier_wrapper(T,train,test,features,pca_perc,...
                            is,just_spikes,rm_non_temporal,which_refs{ir});
            
                    end
        
                    % Do the perfcurve
                    [out.X,out.Y,~,out.AUC] = perfcurve(out.class,out.scores,out.pos_class);
        
                    % Fill the appropriate structure entry
                    model(im).val(iv).side(is).result = out;
        
                end
                
            end
        end
        
        
        approach(ia).model = model;
        all.approach(ia) = approach(ia);
        all.which_ref = which_refs{ir};
        save([plot_folder,file_name],'all')
        
        
        % test - show results
        if 0
            figure
            tiledlayout(2,3)
            for iv = 1:2
                for im = 1:nmodels
                    nexttile
                    ll = plot(model(im).val(iv).side(1).result.X,...
                        model(im).val(iv).side(1).result.Y,'linewidth',2);
                    hold on
                    lr = plot(model(im).val(iv).side(2).result.X,...
                        model(im).val(iv).side(2).result.Y,'linewidth',2);
                    plot([0 1],[0 1],'k--','linewidth',2)
                    xlabel('False positive rate')
                    ylabel('True positive rate')
                    legend([ll,lr],{sprintf('%s: AUC = %1.2f',...
                        model(im).val(iv).side(1).description,model(im).val(iv).side(1).result.AUC),...
                        sprintf('%s: AUC = %1.2f',...
                        model(im).val(iv).side(2).description,model(im).val(iv).side(2).result.AUC)},'fontsize',15,...
                        'location','southeast')
                    title(sprintf('%s %s',model(im).type,model(im).val(iv).description))
                    set(gca,'fontsize',15)
        
        
                end
            end
        end
    
    end
    
    %% Do subsampling analysis
    fprintf('done, took %1.1f seconds',toc);
    fprintf('\nStarting subsampling analysis.\n')
    tic
    
    % Run the lr_mt to extract features
    [T,features,way,dur,sample,ss,durations] =  lr_mt_multitime([2 3]); 
    empty_class = cellfun(@isempty,T.soz_lats);
    T(empty_class,:) = [];

    % Establish HUP and MUSC as training and testing, respectively
    train  = contains(T.names,'HUP');
    test  = contains(T.names,'MP');
    
    % Restrict to car spikes
    car_spikes = contains(features,sprintf('spikes_%s',which_refs{ir}));
    features = features(car_spikes);
    just_spikes = 1; % Just spikes
    
    % Establish what I am varying
    all_durs = unique(dur);
    ndurs = length(all_durs);
    
    all_samples = unique(sample);
    nsamples = length(all_samples);
    
    all_ss = unique(ss);
    nss = length(all_ss);
    
    % Initialize subsampling data
    cv_data = nan(nss,2,ndurs,nsamples); 
    ext_data = nan(nss,2,ndurs,nsamples); 
    
    % Loop over nss
    for iss = 1:nss
        
        curr_ss = all_ss(iss);
    
        cv_auc_l = nan(ndurs,nsamples);
        cv_auc_r = nan(ndurs,nsamples);
    
        ext_auc_l = nan(ndurs,nsamples);
        ext_auc_r = nan(ndurs,nsamples);
        
        % Loop over ndurs -  these will be different error bar points
        for id = 1:ndurs
            fprintf('\nDoing ss %d, dur %d...',iss,id);
            curr_dur = all_durs(id);
            
            
            % Loop over samples
            for is = 1:nsamples
                curr_sample = all_samples(is);
    
                % Get the relevant features
                relevant_features = contains(features,'_way1_') & ...
                    contains(features,sprintf('_ss%d_',curr_ss)) & ...
                    contains(features,sprintf('_dur%d_',curr_dur)) & ...
                    contains(features,sprintf('_samp%d_',curr_sample));
                curr_features = features(relevant_features);
    
                % Run the internal CV models
                left_int = classifier_wrapper(T(train,:),curr_features,...
                    pca_perc,1,just_spikes,rm_non_temporal,[],which_refs{ir});
                right_int = classifier_wrapper(T(train,:),curr_features,...
                    pca_perc,2,just_spikes,rm_non_temporal,[],which_refs{ir});
    
                % Get ROC stats
                [~,~,~,AUCL] = perfcurve(left_int.class,left_int.scores,left_int.pos_class);
                [~,~,~,AUCR] = perfcurve(right_int.class,right_int.scores,right_int.pos_class);
                cv_auc_l(id,is) = AUCL;
                cv_auc_r(id,is) = AUCR;
    
                % Run the external validation models
                left_ext = validation_classifier_wrapper(T,train,test,curr_features,...
                    pca_perc,1,just_spikes,rm_non_temporal,which_refs{ir});
                right_ext = validation_classifier_wrapper(T,train,test,curr_features,...
                    pca_perc,2,just_spikes,rm_non_temporal,which_refs{ir});
    
                % Get ROC stats
                [~,~,~,AUCL] = perfcurve(left_ext.class,left_ext.scores,left_ext.pos_class);
                [~,~,~,AUCR] = perfcurve(right_ext.class,right_ext.scores,right_ext.pos_class);
                ext_auc_l(id,is) = AUCL;
                ext_auc_r(id,is) = AUCR;
    
            end
            fprintf('median left internal AUC is %1.2f\n',nanmedian(cv_auc_l(id,:)))
    
        end
    
        % save data
        cv_data(iss,1,:,:) = cv_auc_l;
        cv_data(iss,2,:,:) = cv_auc_r;
    
        ext_data(iss,1,:,:) = ext_auc_l;
        ext_data(iss,2,:,:) = ext_auc_r;
    
    end
    fprintf('\nDone with subsampling analysis, took %1.1f seconds.\n',toc)
    
    all.cv_ss = cv_data;
    all.ext_ss = ext_data;
    all.durations = durations;
    
    save([plot_folder,file_name],'all')

end

end