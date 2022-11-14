function prep_data_outcome_model

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
int_folder = [results_folder,'analysis/intermediate/'];
out_folder = [results_folder,'analysis/new_outcome/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% Initialize result variables
all_corrs_sp = cell(npts,1);
all_corrs_pear = cell(npts,1);
all_avg_corrs_sp = nan(npts,1);
all_avg_corrs_pear = nan(npts,1);
all_soz_locs = cell(npts,1);
all_soz_lats = cell(npts,1);
all_soz_bin = cell(npts,1);
all_spikes = cell(npts,1);
all_rl = cell(npts,1);
all_ns = cell(npts,1);
all_labels = cell(npts,1);
all_names = cell(npts,1);
all_fc = cell(npts,1);
all_locs = cell(npts,1);
all_coh = cell(npts,1);
all_bp = cell(npts,1);
all_stereo = nan(npts,1);
all_good_spikes = nan(npts,1);
all_two_year_ilae = cell(npts,1);
all_two_year_engel = cell(npts,1);
all_one_year_ilae = cell(npts,1);
all_one_year_engel = cell(npts,1);
all_surgery = cell(npts,1);
all_anatomy = cell(npts,1);
all_aal = cell(npts,1);
all_brainnetome = cell(npts,1);
all_resec_lat = cell(npts,1);
all_resec_loc = cell(npts,1);
all_ablate_lat = cell(npts,1);
all_ablate_loc = cell(npts,1);
all_bipolar_fc = cell(npts,1);
all_bipolar_locs = cell(npts,1);
all_bipolar_labels = cell(npts,1);
all_bipolar_coh = cell(npts,1);
all_bipolar_bp = cell(npts,1);
all_bipolar_aal = cell(npts,1);
all_bipolar_brainnetome = cell(npts,1);
all_native_locs = cell(npts,1);
all_native_bipolar_locs = cell(npts,1);
%all_stitched_coh = cell(npts,1);
%all_stitched_coh_bi = cell(npts,1);


%% Loop over patients
for p = 1:npts

    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    ns =summ.ns_car;
    spikes = summ.spikes;
    rl = summ.rl;
    labels = summ.labels;
    soz_loc = summ.soz.loc;
    soz_lat = summ.soz.lat;
    soz_labels = summ.soz.labels;
    name = summ.name;
    fc = summ.avg_fc;
    locs = summ.locs;
    coh = wrap_or_unwrap_adjacency_fc_toolbox(summ.avg_coh);
    clinical = summ.clinical;
    all_stereo(p) = clinical.stereo;
    good_spikes = summ.good_spikes;
    all_good_spikes(p) = good_spikes;
    anatomy = summ.anatomy;
    bp = summ.bp;
    bipolar_locs = summ.bipolar_locs;
    bipolar_labels = summ.bipolar_labels;
    bipolar_fc = summ.avg_fc_bi;
    bipolar_coh = summ.avg_coh_bi;
    bipolar_bp = summ.bp_bi;
    native_locs = summ.native_locs;
    native_bipolar_locs = summ.native_locs;
    %stitched_coh = summ.stitched_coh;
    %stitched_coh_bi = summ.stitched_coh_bi;
    
    all_names{p} = name;
    
    %% Make spikes nans if bad spikes
    if ~good_spikes
        spikes = nan(size(spikes));
    end
    
    %% Get soz
    all_soz_locs{p} = soz_loc;
    all_soz_lats{p} = soz_lat;
    
    % binary SOZ identification
    soz_bin = zeros(length(labels),1);
    soz_bin(ismember(labels,soz_labels)) = 1;

    %% Unwrap bipolar coh
    bipolar_coh = wrap_or_unwrap_adjacency_fc_toolbox(bipolar_coh);
    
    %% Remove non-intracranial
    ekg = find_non_intracranial(labels);
    ns(ekg,:) = [];
    spikes(ekg,:) = [];
    labels(ekg) = [];
    soz_bin(ekg) = [];
    rl(ekg,:) = [];
    fc(ekg,:) = []; fc(:,ekg) = [];
    locs(ekg,:) = [];
    coh(ekg,:,:) = [];
    coh(:,ekg,:) = [];
    anatomy(ekg) = [];
    all_anatomy{p} = anatomy;
    bp(ekg,:,:) = [];
    bipolar_locs(ekg,:) = [];
    bipolar_labels(ekg) = [];
    bipolar_fc(ekg,:)=[]; bipolar_fc(:,ekg) = [];
    bipolar_coh(ekg,:,:) = []; bipolar_coh(:,ekg,:) = [];
    bipolar_bp(ekg,:,:) = [];
    native_locs(ekg,:) = [];
    native_bipolar_locs(ekg,:) = [];
    %stitched_coh(ekg,:,:) = [];
    %stitched_coh_bi(ekg,:,:) = [];

    
    %% SOZ bin
    all_soz_bin{p} = soz_bin;
    all_labels{p} = labels;
    
    %% Atlas parcellations
    % AAL
    aal = find_atlas_parcellations(locs,'aal');
    all_aal{p} = aal.enames;
    
    % Brainnetome
    brainnetome = find_atlas_parcellations(locs,'brainnetome');
    all_brainnetome{p} = brainnetome.enames;

    %% Bipolar atlas parcelations
    bi_aal = find_atlas_parcellations(bipolar_locs,'aal');
    all_bipolar_aal{p} = bi_aal.enames;

    bi_brainnetome = find_atlas_parcellations(bipolar_locs,'brainnetome');
    all_bipolar_brainnetome{p} = bi_brainnetome.enames;
    
    %% Compare atlas parcellations to clinical anatomical localizations
    if 0
    table(anatomy,aal.enames,brainnetome.enames)
    end
    
    %% Get averages over times
    avg_spikes = nanmean(spikes,2);
    avg_ns = nanmean(ns,2);
    avg_rl = nanmean(rl,2);
    
    all_spikes{p} = avg_spikes;
    all_ns{p} = avg_ns;
    all_rl{p} = avg_rl;
    all_fc{p} = fc;
    all_locs{p} = locs;
    all_coh{p} = coh;
    all_bp{p} = nanmean(bp,3);
    all_bipolar_fc{p} = bipolar_fc;
    all_bipolar_locs{p} = bipolar_locs;
    all_bipolar_labels{p} = bipolar_labels;
    all_bipolar_bp{p} = nanmean(bipolar_bp,3);
    all_bipolar_coh{p} = bipolar_coh;
    all_native_locs{p} = native_locs;
    all_native_bipolar_locs{p} = native_bipolar_locs;
    %all_stitched_coh{p} = stitched_coh;
    %all_stitched_coh_bi{p} = stitched_coh_bi;
    
    %% Correlate average ns and spikes
    avg_corr_sp = corr(avg_spikes,avg_ns,'rows','pairwise','type','spearman');
    avg_corr_pear = corr(avg_spikes,avg_ns,'rows','pairwise','type','pearson');
    all_avg_corrs_sp(p) = avg_corr_sp;
    all_avg_corrs_pear(p) = avg_corr_pear;
    
    %% Correlate time varying
    ntimes = size(spikes,2);
    corr_sp = nan(ntimes,1);
    corr_pear = nan(ntimes,1);
    for it = 1:ntimes
        corr_sp(it) = corr(spikes(:,it),ns(:,it),'rows','pairwise','type','spearman');
        corr_pear(it) = corr(spikes(:,it),ns(:,it),'rows','pairwise','type','pearson');
    end
    
    all_corrs_sp{p} = corr_sp;
    all_corrs_pear{p} = corr_pear;
    
    %% Outcome and surgery
    all_one_year_ilae{p} = clinical.ilae{1};
    all_one_year_engel{p} = clinical.engel{1};
    all_two_year_ilae{p} = clinical.ilae{2};
    all_two_year_engel{p} = clinical.engel{2};
    all_surgery{p} = clinical.surgery;
    all_resec_lat{p} = clinical.resection_lat;
    all_resec_loc{p} = clinical.resection_loc;
    all_ablate_lat{p} = clinical.ablation_lat;
    all_ablate_loc{p} = clinical.ablation_loc;

end

out.avg_corr_sp = all_avg_corrs_sp;
out.avg_corr_pear = all_avg_corrs_pear;
out.all_corrs_sp = all_corrs_sp;
out.all_corrs_pear = all_corrs_pear;
out.all_soz_locs = all_soz_locs;
out.all_soz_lats = all_soz_lats;
out.all_soz_bin = all_soz_bin;
out.all_spikes = all_spikes;
out.all_locs = all_locs;
out.all_rl = all_rl;
out.all_fc = all_fc;
out.all_ns = all_ns;
out.all_coh = all_coh;
out.all_labels = all_labels;
out.all_names = all_names;
out.all_stereo = all_stereo;
out.good_spikes = all_good_spikes;
out.all_one_year_ilae = all_one_year_ilae;
out.all_one_year_engel = all_one_year_engel;
out.all_two_year_ilae = all_two_year_ilae;
out.all_two_year_engel = all_two_year_engel;
out.all_surgery = all_surgery;
out.all_resec_lat = all_resec_lat;
out.all_resec_loc = all_resec_loc;
out.all_ablate_lat = all_ablate_lat;
out.all_ablate_loc = all_ablate_loc;
out.all_anatomy = all_anatomy;
out.all_aal = all_aal;
out.all_brainnetome = all_brainnetome;
out.aal_names = aal.atlas_names;
out.brainnetome_names = brainnetome.atlas_names;
out.all_bp = all_bp;
out.all_bipolar_fc = all_bipolar_fc;
out.all_bipolar_locs = all_bipolar_locs;
out.all_bipolar_bp = all_bipolar_bp;
out.all_bipolar_coh = all_bipolar_coh;
out.all_bipolar_aal = all_bipolar_aal;
out.all_bipolar_brainnetome = all_bipolar_brainnetome;
out.all_native_bipolar_locs = all_native_bipolar_locs;
out.all_native_locs = all_native_locs;
out.all_bipolar_labels = all_bipolar_labels;
%out.all_stitched_coh = all_stitched_coh;
%out.all_stitched_coh_bi = all_stitched_coh_bi;


%% Save
save([out_folder,'main_out.mat'],'out');

end