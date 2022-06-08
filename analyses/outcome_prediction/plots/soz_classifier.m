function soz_classifier(which_atlas,nb)

%{
This function models the SOZ based on various features
%}

%% Initial parameters
% nb = How many random testing/training splits (do 1,000 for paper)
do_plot = 0;
params.models = {'ana','ana_cov','ana_cov_ns','ana_cov_spikes','ana_cov_spikes_ns',}; % which models
params.pretty_name = {'Anatomy','Add coverage density','Add connectivity','Add spike rates','All'};


%% File locations
locations = fc_toolbox_locs;
nmodels = length(params.models);
all_auc = nan(nb,nmodels);

% initialize fun model facts
for im = 1:nmodels
    model_info(im).name = params.models{im};
    model_info(im).pretty_name = params.pretty_name{im};
    model_info(im).all_auc = nan(nb,1);
    model_info(im).all_roc = cell(nb,1);
end

%% Load density model (this is how I get the search radius)
mout = load([locations.paper_plot_folder,'dens_model.mat']);
mout = mout.out;
resid = mout.resid;
params.resid = resid;


params.which_atlas = which_atlas;
params.sr = mout.sr; % search radius (leave empty to use default calc)
params.prop_train = 2/3;
params.do_r2 = 0; % r^2 instead of r for FC measurement? (should be zero)
params.do_ns_resid = 0; % take residuals of ns (density normalized).
params.include_lat = 0; % include laterality in addition to broad anatomical regions. I don't do this but should I????
params.dens_model = 1; % use Erin's density model as opposed to rat11. Doesn't matter because I don't do residuals.

%params.max_spikes = 1/3600; % max spikes for both distance and density models
params.do_norm = 1; % normalize spike rate and density within pt
params.atlas_anatomy = 0; % break anatomy into many categories (all atlas regions)? Model fails to converge
%params.only_sleep = 0; % outside scope, don't change
params.test_implant = 0; % 0: all patients, 1: sEEG, 2: grid/strip


%% Get AUC for each model for each random testing/training split
params.do_glme = 0; 
for ib = 1:nb
    if mod(ib,10) == 0, fprintf('\nDoing %d of %d\n',ib,nb); end
    out = individual_classifier(params);
    for im = 1:nmodels
        all_auc(ib,im) = out.model(im).auc;
        model_info(im).all_auc(ib) = out.model(im).auc;
        model_info(im).all_roc{ib} = out.model(im).roc;
    end
end
all_out.npts_remain = out.npts_remain;

%% for each pair of models, generate stats for AUC difference 
model_p = nan(nmodels,nmodels);
for im = 1:nmodels
    for jm = 1:im-1
        
        % Get aucs of each
        auc_i = model_info(im).all_auc;
        auc_j = model_info(jm).all_auc;
        
        % bootstrap stats
        stats = bootstrap_ci_and_p(auc_i,auc_j);
        model_p(im,jm) = stats.p;
        
    end
end

%% Also estimate bootstap CI of own AUC
bootci = nan(nmodels,2);
for im = 1:nmodels
    auc = model_info(im).all_auc;
    bootci(im,:) = [prctile(auc,2.5) prctile(auc,97.5)];
end

%array2table(all_auc,'variablenames',models)
for im = 1:nmodels
    [unify_x,unify_y] = unify_roc(model_info(im).all_roc);
    model_info(im).x = unify_x;
    model_info(im).y = unify_y;
    model_info(im).ym = mean(unify_y);
    model_info(im).ly = prctile(unify_y,5,1);
    model_info(im).uy = prctile(unify_y,95,1);
    
end

all_out.model_info = model_info;
all_out.model_p = model_p;
all_out.bootci = bootci;

if do_plot
    figure
    set(gcf,'position',[10 10 1400 350])
    tiledlayout(1,5,'tilespacing','tight','padding','tight')

    for im = [1 2 3 4 5]
        nexttile
        mp = shaded_error_bars_fc(model_info(im).x,model_info(im).ym,...
            [model_info(im).ly',model_info(im).uy'],'k');
        hold on
        plot([0 1],[0 1],'k--','linewidth',2)
        title(model_info(im).pretty_name)
        legend(mp,sprintf('Mean AUC = %1.2f',median( model_info(im).all_auc)),...
            'location','southeast','fontsize',15)
        set(gca,'fontsize',15)
        xlabel('False positive rate')
        ylabel('True positive rate')
    end

    print(gcf,[locations.paper_plot_folder,'prediction_nospikes'],'-dpng')
end


%% Do once with all data to get glme paramaters
params.do_glme = 1; 
glme_beta = nan(nb,1);
for ib = 1:nb
    glme_stuff = individual_classifier(params);
    beta = glme_stuff.model(5).glm.Coefficients{5,2};
    glme_beta(ib) = beta;
end
all_out.glme_stuff.last_model = out;
all_out.glme_stuff.all_betas = glme_beta;
% Get stats on the betas
stats = bootstrap_ci_and_p(glme_beta);
all_out.glme_stuff.stats = stats;


%% Do it again with stereo vs grid strip
stereo_all_auc = nan(nb,nmodels);
not_stereo_all_auc = nan(nb,nmodels);
params.do_glme = 0; % reset to training/testing paradigm
for ib = 1:nb
    if mod(ib,10) == 0, fprintf('\nDoing %d of %d\n',ib,nb); end
    params.test_implant = 1;
    stereo_out = individual_classifier(params);
    
    params.test_implant = 2;
    not_stereo_out = individual_classifier(params);
    for im = 1:nmodels
        stereo_all_auc(ib,im) = stereo_out.model(im).auc;
        stereo.model_info(im).all_auc(ib) = stereo_out.model(im).auc;
        stereo.model_info(im).all_roc{ib} = stereo_out.model(im).roc;
        
        not_stereo_all_auc(ib,im) = not_stereo_out.model(im).auc;
        not_stereo.model_info(im).all_auc(ib) = not_stereo_out.model(im).auc;
        not_stereo.model_info(im).all_roc{ib} = not_stereo_out.model(im).roc;
    end
end
stereo_vs_not.n_stereo_remain = stereo_out.n_stereo;
stereo_vs_not.n_not_stereo_remain = not_stereo_out.n_not_stereo;

% for each model, compare stereo to not stereo
for im = 1:nmodels
    
    % Get aucs of each
    stereo_auc = stereo.model_info(im).all_auc;
    not_stereo_auc = not_stereo.model_info(im).all_auc;

    % Get stats on the diff
    stats = bootstrap_ci_and_p(stereo_auc,not_stereo_auc);
        
    stereo_vs_not.model(im).p = stats.p;
    stereo_vs_not.model(im).means = [nanmean(stereo_auc),nanmean(not_stereo_auc)];
    stereo_vs_not.model(im).ci95 = [stats.CI_95];
    stereo_vs_not.model(im).all = [stereo_auc',not_stereo_auc'];
  
end

all_out.stereo_vs_not = stereo_vs_not;

save([locations.paper_plot_folder,'model_stuff_',params.which_atlas,'.mat'],'all_out')



end

function mout = individual_classifier(params)
%{
Predict whether an electrode is in the SOZ based on 1) its spike rate
normalized across other electrodes within the patient and 2) its anatomical
localization (a categorical variable that can be MT, TN, OC, or nothing)
%}



%% Parameters
%max_spikes = params.max_spikes; % max spikes/elecs/min to include in model
which_atlas = params.which_atlas;
prop_train = params.prop_train;
do_norm = params.do_norm;
models = params.models;
atlas_anatomy = params.atlas_anatomy;
do_ns_resid = params.do_ns_resid;
do_glme = params.do_glme;
include_lat = params.include_lat;
dens_model = params.dens_model;
do_r2 = params.do_r2;
sr = params.sr;
test_implant = params.test_implant;

%% File locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
spike_folder = locations.paper_data_folder;
atlas_folder = locations.paper_data_folder;


%% Load out file and get roc stuff
out = load([spike_folder,'main_out.mat']);
out = out.out;

%% Get stuff
rate = out.all_spikes;
soz = out.all_soz_bin;
npts = length(soz);
labels = out.all_labels;
fc = out.all_fc;
locs = out.all_locs;
all_stereo = out.all_stereo;
stereo = all_stereo == 1; n_stereo = sum(stereo);
not_stereo = all_stereo == 0; n_not_stereo = sum(not_stereo);

%% Turn soz to logical
soz = cellfun(@logical,soz,'uniformoutput',false);

%% R^2???
if do_r2
    fc = cellfun(@(x) x.^2, fc,'uniformoutput',false);
end

%% Load atlas file
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;

%% get atlas stuff
atlas_elec_labels = atlas_out.elecs_labels;
atlas_elec_regions = atlas_out.elecs_atlas;
atlas_nums = atlas_out.atlas_nums;
atlas_names = atlas_out.atlas_names;

%% Spatially normalize the FC matrix
%
if dens_model
    resid = params.resid;
else
    resid = fit_distance_model(locs,fc,soz,rate,max_spikes,plot_folder);
end

ns_resid = cellfun(@(x) nansum(x,2),resid,'uniformoutput',false);
ns = cellfun(@(x) nanmean(x,2),fc,'uniformoutput',false); % average rather than sum

if do_ns_resid == 1
    ns = ns_resid;
end




%% Localize regions into broad categories
broad = localize_regions(atlas_names,which_atlas);
if ~include_lat
    broad_no_lat = cell(length(broad),1);
    for ib = 1:length(broad)
        curr = broad{ib};
        if isempty(curr),continue; end
        curr = strrep(curr,'left ','');
        curr = strrep(curr,'right ','');
        broad_no_lat{ib} = curr;
    end
else
    broad_no_lat = broad;
end

%% Get localizations for each patient
elec_broad = cell(npts,1);
for ip = 1:npts
    elec_broad{ip} = get_anatomy_elecs(labels{ip},atlas_elec_labels{ip},atlas_elec_regions{ip},broad_no_lat,atlas_nums,atlas_anatomy);
end

if 0
    table(elec_broad{110},labels{110})
end

%% Get stuff into single column format
vec_rate = [];
vec_loc = {};
vec_soz = [];
vec_pt_idx = [];
vec_ns = [];
vec_dens = [];

%% GEt default search radius
if isempty(sr)
    sr = calculate_default_search_radius(locs); 
    %sr = max_inter_elec_dist(locs);
end

for ip = 1:npts
    curr_rate = rate{ip};
    curr_soz = soz{ip};
    curr_ana_loc = elec_broad{ip};
    curr_ns = ns{ip};
    curr_loc = locs{ip};
    
    if atlas_anatomy
        curr_ana_loc = arrayfun(@(x) num2str(x),curr_ana_loc,'uniformoutput',false);
    end
    
    % get density
    density = estimate_coverage_density(curr_loc,sr);
        
    % normalize within patient
    if do_norm
        curr_rate = (curr_rate-nanmean(curr_rate))./nanstd(curr_rate);
        %
        if ~do_ns_resid
            curr_ns = (curr_ns-nanmean(curr_ns))./nanstd(curr_ns);
        end
        
        density = (density - nanmean(density))./nanstd(density);
    end
    
    vec_rate = [vec_rate;curr_rate];
    vec_dens = [vec_dens;density];
    vec_loc = [vec_loc;curr_ana_loc];
    vec_soz = [vec_soz;curr_soz];
    vec_ns = [vec_ns;curr_ns];
    vec_pt_idx = [vec_pt_idx;repmat(ip,length(curr_soz),1)];
end

%% Make table
T = table(vec_soz,vec_pt_idx,vec_rate,vec_loc,vec_ns,vec_dens);


%% Remove rows without localization
assert(length(unique(T.vec_pt_idx))==npts); % assure I start with all patients
no_loc = cellfun(@isempty,T.vec_loc);

T(no_loc,:) = [];

% number excluded due to localization
n_exc_locs = npts - length(unique(T.vec_pt_idx)); % subtract those who remain to find those excluded


%% Remove nan rows % patients with bad spikes or other weirdness
nan_rows = isnan(T.vec_soz) | isnan(T.vec_pt_idx) | isnan(T.vec_rate) | ...
    isnan(T.vec_ns) | isnan(T.vec_dens);

% confirm only excluding due to spikes
assert(sum(isnan(T.vec_soz) | isnan(T.vec_pt_idx)|isnan(T.vec_ns) | isnan(T.vec_dens)) == 0);
T(nan_rows,:) = [];

% number excluded due to spikes
n_exc_spikes = npts - n_exc_locs - length(unique(T.vec_pt_idx));


%% Get number of patients
unique_pts = unique(T.vec_pt_idx);
npts_remain = length(unique_pts);
assert(npts == npts_remain+n_exc_spikes+n_exc_locs);
mout.npts = npts;
mout.npts_remain = npts_remain;
mout.n_exc_locs = n_exc_locs;
mout.n_exc_spikes = n_exc_spikes;

%% Divide into training and testing set
while 1 % wrap this all in a while loop to try again if wacky errors related to bad selection (occurs in ~1/100-1/200 times)
    continue_out_of_while_loop = 0;
    if do_glme
        npts = length(unique_pts);
        training = randsample(unique_pts,npts,true); % the true means to sample WITH replacement

        % rebuild the table
        temp_table  = cell2table(cell(0,size(T,2)),'VariableNames',T.Properties.VariableNames);
        for ip = 1:length(training)
            % get the rows of T with this patient id
            curr_rows = ismember(T.vec_pt_idx,training(ip));

            temp_table = [temp_table;T(curr_rows,:)];

        end

        T_train = temp_table;
        T_test = temp_table;

        assert(isequal(unique(T_train.vec_pt_idx),unique(training)))
        %training_idx = ismember(T.vec_pt_idx,training);
        %T_train = T(training_idx,:);
    else

        if test_implant == 1 % stereo

            stereo_idx = stereo(T.vec_pt_idx) == 1; % Reduce to only stereo
            T(~stereo_idx,:) = [];
            unique_pts = unique(T.vec_pt_idx);
            mout.n_stereo = length(unique_pts);
            training = randsample(n_stereo,floor(n_stereo*prop_train));
            training_idx = ismember(T.vec_pt_idx,training);
            testing_idx = ~ismember(T.vec_pt_idx,training);
            T_train = T(training_idx,:);
            T_test = T(testing_idx,:);

            % confirm all stereo
            assert(~any(stereo(T.vec_pt_idx) ~= 1));

        elseif test_implant == 2 % grid/strip

            not_stereo_idx = not_stereo(T.vec_pt_idx) == 1; % Reduce to only not stereo
            T(~not_stereo_idx,:) = [];
            unique_pts = unique(T.vec_pt_idx);
            mout.n_not_stereo = length(unique_pts);
            training = randsample(n_not_stereo,floor(n_not_stereo*prop_train));
            training_idx = ismember(T.vec_pt_idx,training);
            testing_idx = ~ismember(T.vec_pt_idx,training);
            T_train = T(training_idx,:);
            T_test = T(testing_idx,:);

            assert(~any(not_stereo(T.vec_pt_idx) ~= 1));

        else

            training = randsample(npts,floor(npts*prop_train));
            training_idx = ismember(T.vec_pt_idx,training);
            testing_idx = ~ismember(T.vec_pt_idx,training);
            T_train = T(training_idx,:);
            T_test = T(testing_idx,:);


        end

        %% Confirm that I separated by patients
        train_pts = unique(T_train.vec_pt_idx);
        test_pts = unique(T_test.vec_pt_idx);
        assert(isempty(intersect(train_pts,test_pts)))

    end



    % Loop over models
    for im = 1:length(models)
        curr_model = models{im};
        switch curr_model

            case 'ana'
                formula = 'vec_soz ~ vec_loc';
                formula_me = 'vec_soz ~ vec_loc + (1|vec_pt_idx)';
            case 'ana_cov'
                formula = 'vec_soz ~ vec_loc + vec_dens';
                formula_me = 'vec_soz ~ vec_loc + vec_dens+ (1|vec_pt_idx)';
            case 'ana_cov_spikes'
                formula = 'vec_soz ~ vec_loc + vec_dens + vec_rate';
                formula_me = 'vec_soz ~ vec_loc + vec_dens + vec_rate + (1|vec_pt_idx)';
            case 'ana_cov_spikes_ns'
                formula = 'vec_soz ~ vec_loc + vec_dens + vec_rate + vec_ns';
                formula_me = 'vec_soz ~ vec_loc + vec_dens + vec_rate + vec_ns + (1|vec_pt_idx)';
            case 'ana_cov_ns'
                formula = 'vec_soz ~ vec_loc + vec_dens + vec_ns';
                formula_me = 'vec_soz ~ vec_loc + vec_dens + vec_ns + (1|vec_pt_idx)';
        end

        %% Do both the glm and glme
        %T_train_me = T_train;
       % T_train_me.vec_pt_idx = nominal(T_train_me.vec_pt_idx);
       
       try
           if do_glme % try glme regardless
               T_train.vec_pt_idx = nominal(T_train.vec_pt_idx);
               T_test.vec_pt_idx = nominal(T_test.vec_pt_idx);
               glm = fitglme(T_train,formula_me,'Distribution','Binomial');
           else
               glm = fitglm(T_train,formula,'Distribution','Binomial');
           end
       catch ME
            if contains(ME.message,'NaN or Inf values are not allowed in X')
                fprintf('\nWacky error, resample again...\n');
                continue_out_of_while_loop = 1;
                break % break out of model for loop
            else
                error('unexpected error');
            end
           
       end
      %  

        %% Test glm on testing data
        params = glm.CoefficientNames;

        if do_glme
            mout.model(im).params = params;
            mout.model(im).name = curr_model;
            mout.model(im).glm = glm;
            continue;
        end

        if 0
            asum = zeros(size(T_test,1),1);
            asum = asum + glm.Coefficients.Estimate(1);
            for p = 2:length(params)
                est = glm.Coefficients.Estimate(p);
                if contains((params{p}),'vec_loc')
                    C = strsplit((params{p}),'_');
                    cat = C{3};

                    % find the indices of T_test that fit this category
                    match = strcmp(T_test.vec_loc,cat);
                    asum(match) = asum(match) + est;
                else

                    asum = asum +  T_test.(params{p})*est;
                end
            end


            classification = logistic(asum);
        else
            % Derive classifications (probability of SOZ) for testing data
            %classification = feval(glm,T_test); % predict
            classification = predict(glm,T_test);
            
        end

        % Compare classification against true
        [X,Y,~,AUC,~] = perfcurve(T_test.vec_soz,classification,1);

        mout.model(im).name = curr_model;
       % out.model(im).glme = glme;
        mout.model(im).roc = [X,Y];
        mout.model(im).auc = AUC;

        if isempty(mout.model(im).auc)
            error('what');
        end
    end
    if continue_out_of_while_loop == 1
        continue;
    end
    break % if survive to here, break out of the while loop
end





end

function elec_broad = get_anatomy_elecs(sp_labels,atlas_labels,regions,broad_no_lat,atlas_nums,atlas_anatomy)
    
%% Remove ekg from atlas
ekg = find_non_intracranial(atlas_labels);
atlas_labels(ekg) = [];
regions(ekg) = [];

%% Reconcile labels
assert(isequal(sp_labels,atlas_labels)) 

if atlas_anatomy
    elec_broad = regions;
    return
end

elec_broad = cell(length(sp_labels),1);

% Loop over electrodes
for ich = 1:length(sp_labels)
    
    % get the region of this electrode
    curr_region = regions(ich);
    
    % get the index of this
    curr_idx = find(atlas_nums == curr_region);
    
    if isempty(curr_idx), continue; end
    
    elec_broad{ich} = broad_no_lat{curr_idx};
end

if 0
    table(sp_labels,elec_broad)
end


end

function out = logistic(x)

out = 1./(1+exp(-x));

end


