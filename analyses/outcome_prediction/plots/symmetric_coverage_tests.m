function nout = symmetric_coverage_tests(which_atlas,force_same_num)

%{
This performs the analysis comparing connectivity to the SOZ to the region
contralateral to the SOZ, first restricting regions to only those with
bilateral coverage.

Double checked 6/22/22, looks good.
%}

%% Parameters
no_symm_cov = 0; % set to 1 if you want to see the results if you DON'T restrict to symmetric coverage (should be 0 for paper)
randomize_lat = 0; % make fake SOZ and SOZ lats (as a test to confirm negative result if I randomize lats)
do_plots = 0; 
freqs = {'Delta (0.5-4 Hz)','Theta (4-8 Hz)','Alpha (8-12 Hz)','Beta (12-30 Hz)','Gamma (30-80 Hz)'};

%% Get file locs
locations = fc_toolbox_locs;
plot_folder = locations.paper_plot_folder;

atlas_folder = locations.paper_data_folder;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load out file with functional connectivity and spikes as well as SOZ info
out = load([locations.paper_data_folder,'main_out.mat']);
out = out.out;

%% Get stuff
nfreqs = length(freqs);
rate = out.all_spikes;
soz = out.all_soz_bin;
npts = length(soz);
labels = out.all_labels;
fc = out.all_fc;
coh = out.all_coh;

%% Get total n for main result n
any_locs = zeros(npts,1);
for ip = 1:npts
   if any(~isnan(out.all_locs{ip}),'all')
       any_locs(ip) = 1;
   end
end
nout.pts_with_any_locs = (any_locs);

% Soz lats
soz_lats = out.all_soz_lats;
right_lat = strcmp(soz_lats,'right');
left_lat = strcmp(soz_lats,'left');

% randomize lats to confirm negative result if randomized
if randomize_lat
    nright = sum(right_lat);
    nleft = sum(left_lat);
    ntotal = length(right_lat);
    
    temp_right = randsample(ntotal,nright);
    right_lat = zeros(ntotal,1);
    right_lat(temp_right) = 1;
    right_lat = logical(right_lat);
    
    temp_left = 1:ntotal;
    temp_left(right_lat) = [];
    temp_left = randsample(temp_left,nleft);
    left_lat = zeros(ntotal,1);
    left_lat(temp_left) = 1;
    left_lat = logical(left_lat);
    
end

%% Load atlas file
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;

assert(isequal(out.all_names,atlas_out.pt_names))

%% get atlas stuff (indictates which electrodes are in which atlas regions)
atlas_elec_labels = atlas_out.elecs_labels;
atlas_elec_regions = atlas_out.elecs_atlas;
atlas_nums = atlas_out.atlas_nums;
names = atlas_out.atlas_names;

%% Get locs and lats for atlas names
[locs,lats] = lateralize_regions(names,which_atlas);
left = strcmp(lats,'L');
right = strcmp(lats,'R');
neither_lat = ~left & ~right;

% confirm atlas has as many on right as left
assert(sum(left)==sum(right));

%% Construct atlas (regions x regions x patients)
% orig_atlas
%{
orig_atlas = rebuild_atlas(fc,rate,atlas_elec_labels,...
    atlas_elec_regions,atlas_nums,labels,soz,coh,0,locs,lats);
%}

% Puts FC, spikes, etc. from electrode space into atlas space
[atlas,spikes,bin_soz,coh] = rebuild_atlas(fc,rate,atlas_elec_labels,...
    atlas_elec_regions,atlas_nums,labels,soz,coh,force_same_num,locs,lats);



%% Compare some elec labels to atlas labels
% making sure that LA is left amygdala, etc...
% did some spot checks, looks ok
if 0
    ip = 1;
    non_nan = ~isnan(atlas_elec_regions{ip});
    atlas_names = names(atlas_elec_regions{ip}(non_nan));
    atlas_labels = atlas_elec_labels{ip}(non_nan);
    table(atlas_labels,atlas_names)
    
end

%% Re-order atlas to be left then right then neither
lr_order = reorder_lr(locs,lats); % get the order

% Put everything into this order
left = left(lr_order);
right = right(lr_order);
neither_lat = neither_lat(lr_order);
atlas = atlas(lr_order,lr_order,:);
%orig_atlas = orig_atlas(lr_order,lr_order,:);
coh = coh(lr_order,lr_order,:,:);
names = names(lr_order);
locs = locs(lr_order);
lats = lats(lr_order);
spikes = spikes(lr_order,:);
bin_soz = bin_soz(lr_order,:);
atlas_nums = atlas_nums(lr_order);

%% index of contralateral region for each region
contra_index = nan(size(left));
contra_index(1:sum(left)) = ([1:sum(left)])'+sum(left); % the right sided ones start after left, so contralateral ones for the first nleft are these
contra_index(sum(left)+1:sum(left)*2) = ([1:sum(left)])';

% make sure first half locs are same as last half locs
assert(isequal(locs(1:sum(left)),locs(sum(left)+1:sum(left)*2)))
assert(isequal(locs(contra_index(1:sum(left))),locs(contra_index(sum(left)+1:sum(left)*2))))

%% Double check contralateral is contralateral
if 0
    % looks good
    nan_contra = isnan(contra_index);
    names_no_contra = names(~nan_contra);
    table(names_no_contra,names_no_contra(contra_index(~nan_contra)))
end

%{
%% Define names corresponding to mesial temporal (not used in paper)
switch which_atlas
    case 'aal_bernabei'
        mt_names = {'Hippocampus','Amygdala'};
    case 'brainnetome'
        mt_names = {'Amyg','Hipp'};
end

mt = contains(names,mt_names);
%}

%% First, build symmetric coverage atlas
% Build atlas for correlation and get indices of bilateral coverage regions
[symm_cov_atlas,all_bilateral] = build_symmetric_coverage_atlas(atlas,locs,lats);
%[symm_cov_atlas_orig,all_bilateral_orig] = build_symmetric_coverage_atlas(orig_atlas,locs,lats);

if ~no_symm_cov
    atlas = symm_cov_atlas; % redefine the FC atlas to be the symmetric coverage atlas
else
    % stuff for the test of allowing non symmetric coverage
    all_bilateral = zeros(size(all_bilateral));
    for ip = 1:npts
        curr_atlas = atlas(:,:,ip);
        avg_row = nanmean(curr_atlas,2);
        all_bilateral(:,ip) = ~isnan(avg_row);
        
    end
end

% Build symmetric coverage coherence matrix
symm_cov_coh = nan(size(coh));
for ip = 1:npts
    curr_bilateral = logical(all_bilateral(:,ip));
    curr_coh = coh(:,:,:,ip);
    curr_coh(~curr_bilateral,:,:) = nan; % set to nan if not bilateral coverage
    curr_coh(:,~curr_bilateral,:) = nan;
    symm_cov_coh(:,:,:,ip) = curr_coh;
end
coh = symm_cov_coh;

% Build symmetric coverage spike matrix
% I believe this is also important for spikes because some regions (like
% hippocampus) have more spikes and are preferentially implanted on the
% contralateral hemisphere as sentinel electrodes. Without restricting to
% symmetric coverage, we would have a bias toward finding MORE spikes on
% the opposite side of the SOZ because this side is mostly hippocampal
% coverage whereas the SOZ side would average a broad less-spikey area.
symm_spikes = nan(size(spikes));
for ip = 1:npts
    curr_bilateral = logical(all_bilateral(:,ip));
    curr_spikes = spikes(:,ip);
    curr_spikes(~curr_bilateral) = nan;
    symm_spikes(:,ip) = curr_spikes;
end
spikes = symm_spikes;

%% Double check symmetric coverage
if ~no_symm_cov
% Looks good!
for ip = 1:npts
    curr = atlas(:,:,ip);
    avg_rows = nanmean(curr,2);
    
    % find non nan
    non_nan = find(~isnan(avg_rows)); % find rows for that patient with any network edges
    
    % confirm contralateral is non nan (make sure contralateral region also
    % has network edges)
    for i = 1:length(non_nan)
        assert(~isnan(avg_rows(contra_index(non_nan(i))))|length(non_nan)==3) % latter is an edge case, forgot why
    end
    
end
end

%% Get number of patients with any symmetric coverage and average number of symmetric regions
nsymmetric = sum(all_bilateral,1)/2; % number (unilateral) for each patient
any_symmetric = sum(nsymmetric>0);
total_n_for_symmetric = length(nsymmetric);
mean_symmetric = mean(nsymmetric);
nout.n_coverage.nsymmetric = nsymmetric;
nout.n_coverage.any_symmetric = any_symmetric;
nout.total_n_for_symmetric = total_n_for_symmetric;
nout.mean_symmetric = mean_symmetric;


%% Get average connectivity of SOZ (to other things) and that of contralateral region
% INitialize things
soz_intra = nan(npts,2); % soz to itself vs contralateral
lr_soz_intra = nan(npts,2); % left soz loc to itself vs right
soz_all = nan(npts,2); % avg connectivity to soz vs contralateral
lr_soz_all = nan(npts,2); % same but L-R
hemi = nan(npts,2); % intrahemispheric connectivity on side of SOZ vs contra
hemi_lr = nan(npts,2); % same but LR
soz_coh_all = nan(npts,nfreqs,2); % avg coherence to SOZ
all_bin_contra_soz = zeros(size(bin_soz)); 
ns_soz_minus_contra = nan(sum(strcmp(lats,'R')),npts);

%% Show common soz locations
if 0 
    npts_soz = sum(bin_soz,2);
    [~,sorted_soz] = sort(npts_soz,'descend')
    table(names(sorted_soz(1:10)))
    % Doing this shows that the most common SOZ regions are left
    % hippocampus and then right hippocampus, which is not surprising
end

for ip = 1:npts
    
    %% SOZ connectivity
    % get soz indices
    curr_soz = bin_soz(:,ip);
    
    if randomize_lat % test
        nsoz = sum(curr_soz);
        k = randsample(length(curr_soz),nsoz);
        curr_soz = zeros(length(curr_soz),1);
        curr_soz(k) = 1;
        curr_soz = logical(curr_soz);
    end
    
    % get the regions contalateral to soz
    contra_soz = contra_index(curr_soz);
    contra_soz(isnan(contra_soz)) = [];
    bin_contra_soz = zeros(length(curr_soz),1);
    bin_contra_soz(contra_soz) = 1;
    bin_contra_soz = logical(bin_contra_soz);
    all_bin_contra_soz(:,ip) = bin_contra_soz;
    
    % Show some
    if 0
        table(names(curr_soz),names(bin_contra_soz))
    end
    
    % confirm that difference in number of soz and contra soz is just those
    % regions that have no contralateral thing. Note, some patients will
    % have bilateral SOZ.
    switch which_atlas
        case 'aal_bernabei' % some regions without laterality
            assert(abs(sum(curr_soz)-sum(bin_contra_soz)) == ...
                sum(ismember([names(curr_soz)';names(bin_contra_soz)'],names(end-8:end))))
        case 'brainnetome'
            assert((sum(curr_soz)==sum(bin_contra_soz)))
    end
    
    
    %% Get the left and right
    all_ipsi_and_contra_soz = [find(curr_soz);find(bin_contra_soz)];
    left_soz = all_ipsi_and_contra_soz(left(all_ipsi_and_contra_soz) == 1);
    right_soz = all_ipsi_and_contra_soz(right(all_ipsi_and_contra_soz) == 1);
    
    assert(sum(strcmp(lats(left_soz),'L')) == length(left_soz) && ...
        sum(strcmp(lats(right_soz),'R'))==length(right_soz));
    
    lr_soz_intra(ip,:) =  [nanmean(atlas(left_soz,left_soz,ip),'all'),...
        nanmean(atlas(right_soz,right_soz,ip),'all')];
    lr_soz_all(ip,:) =  [nanmean(atlas(left_soz,:,ip),'all'),...
        nanmean(atlas(right_soz,:,ip),'all')];
    
    %% get SOZ-SOZ connectivity and contra-SOZ - contra-SOZ connectivity
    soz_intra(ip,:) = [nanmean(atlas(curr_soz,curr_soz,ip),'all'),...
        nanmean(atlas(bin_contra_soz,bin_contra_soz,ip),'all')];
        
    soz_all(ip,:) = [nanmean(atlas(curr_soz,:,ip),'all'),...
        nanmean(atlas(bin_contra_soz,:,ip),'all')];
    
    soz_coh_all(ip,:,1) = nanmean(coh(curr_soz,:,:,ip),[1 2]);
    soz_coh_all(ip,:,2) = nanmean(coh(bin_contra_soz,:,:,ip),[1 2]);
    
    %% Holohemispjheric
    % get everything on the side of the soz
    % Note that here my definition of SOZ is different - it is the final
    % anatomical localization at the end of the EMU discharge summary
    % rather than the seizure-by-seizure individual electrode
    if right_lat(ip) == 1
        ipsi_lats = strcmp(lats,'R');
        contra_lats = strcmp(lats,'L');
    elseif left_lat(ip) == 1
        ipsi_lats = strcmp(lats,'L');
        contra_lats = strcmp(lats,'R');
    else
        continue;
    end
    
    % Get soz side - soz side connectivity and contralateral connectivity
    hemi(ip,:) = [nanmean(atlas(ipsi_lats,ipsi_lats,ip),'all'),...
        nanmean(atlas(contra_lats,contra_lats,ip),'all')];
   
    hemi_lr(ip,:) = [nanmean(atlas(left,left,ip),'all'),...
        nanmean(atlas(right,right,ip),'all')];
    
    %% NS
    ns = nanmean(atlas(:,:,ip),2);
    ipsi = find(ipsi_lats);
    contra = contra_index(ipsi);
    ns_soz_minus_contra(:,ip) = ns(ipsi)-ns(contra);
    
    if 0
        table(names(ipsi)',names(contra)')
    end
end

%% For each region, get % where ns_soz_minus_contra is negative
perc_neg = nan(size(ns_soz_minus_contra,1),1);
for i = 1:length(perc_neg)
    perc_neg(i) = sum(ns_soz_minus_contra(i,:)<0)/sum(~isnan(ns_soz_minus_contra(i,:)));
    
end

if 0
    table((names(left))',perc_neg)
end

if 0
    turn_nans_gray(ns_soz_minus_contra)
    yticks(1:size(ns_soz_minus_contra,1))
    yticklabels(names(ipsi))
    caxis([min(ns_soz_minus_contra,[],'all') -min(ns_soz_minus_contra,[],'all')])
end

if ~no_symm_cov && ~randomize_lat
%% Get numbers for above analyses
n_soz_all = sum(~any(isnan(soz_all),2)); % patients with symmetric coverage of SOZ and contralateral region
nout.n_analyses.n_soz_all = n_soz_all;
% Is above number what I expect it to be?
symm_cov_soz = zeros(npts,1);
for ip = 1:npts
    curr_soz = bin_soz(:,ip);
    curr_atlas = atlas(:,:,ip);
    if sum(~isnan(curr_atlas(curr_soz,:)),'all') ~= 0
        symm_cov_soz(ip) = 1;
    end
end
assert(sum(symm_cov_soz) == n_soz_all)

% Patients with unilateral epilepsy and at least two symmetric regions in
% each hemisphere for intra-hemispheric analyses
n_holo_hem = sum(~any(isnan(hemi),2));
nout.n_analyses.n_holo_hem = n_holo_hem;
% confirm above number is what i expect it to be
unilat_and_2_regions = zeros(npts,1);
for ip = 1:npts
    unilat = right_lat(ip) == 1 || left_lat(ip) == 1;
    avg_row = nanmean(atlas(:,:,ip),2);
    two_regions = sum(~isnan(avg_row)) > 3; % 4 total (2 on each side)
    
    if unilat && two_regions
        unilat_and_2_regions(ip) = 1;
    end
    
end
assert(n_holo_hem == sum(unilat_and_2_regions))

% Patients with at least 2 SOZ regions and contralateral coverage
n_soz_intra = sum(~any(isnan(soz_intra),2));
nout.n_soz_intra = n_soz_intra;
% confirma bove number what I expect it to be
two_soz_and_contra = zeros(npts,1);
for ip = 1:npts
    curr_soz = bin_soz(:,ip);
    avg_row = nanmean(atlas(:,:,ip),2);
    soz_regions = sum(~isnan(avg_row(curr_soz)));
    if sum(soz_regions) >= 2
        two_soz_and_contra(ip) = 1;
    end
    
end
assert(n_soz_intra == sum(two_soz_and_contra));
end

%% Build a SOZ - non SOZ laterality ordered atlas
soz_non_soz_ordered_atlas = build_soz_ordered_atlas(atlas,left,right,right_lat,left_lat);

%% get confusion matrix for connectivity
hemi_diff = hemi_lr(:,1) - hemi_lr(:,2);
predicted = cell(length(soz_lats),1);

% prediction rule: lower connectivity on side of SOZ
predicted(hemi_diff > 0) = {'right'}; % if diff positive, so L more, predict R sided epilepsy
predicted(hemi_diff < 0) = {'left'};
empty = hemi_diff == 0 | isnan(hemi_diff) | (~strcmp(soz_lats,'right') & ~strcmp(soz_lats,'left'));
predicted(empty) = []; % remove those without info
lats_for_conf = soz_lats;
lats_for_conf(empty) = [];

if 0
    table(lats_for_conf,predicted)
end

% Calculate confusion matrix
if force_same_num == 1
    conf_out_fc = [];
else
    conf_out_fc = confusion_matrix(predicted,lats_for_conf,0);
end

% Get numbers for above analysis
n_conf_fc = length(predicted);
if ~randomize_lat
    assert(n_conf_fc == n_holo_hem)
end

%% get confusion matrix for spikes
% Note that this is also using a symmetric coverage map!
hemi_lr_spikes = [nanmean(spikes(left,:),1)',nanmean(spikes(right,:),1)'];
hemi_diff = hemi_lr_spikes(:,1) - hemi_lr_spikes(:,2);
predicted = cell(length(soz_lats),1);

% prediction rule: more spikes on side of SOZ
predicted(hemi_diff < 0) = {'right'};  % if fewer spikes left, predict right sided epilepsy
predicted(hemi_diff > 0) = {'left'};
empty = hemi_diff == 0 | isnan(hemi_diff) | (~strcmp(soz_lats,'right') & ~strcmp(soz_lats,'left'));
predicted(empty) = [];
lats_for_conf = soz_lats;
lats_for_conf(empty) = [];

if 0
    table(lats_for_conf,predicted)
end

if force_same_num == 1
    lat_info = [];
    conf_out_spikes = [];
else
    conf_out_spikes = confusion_matrix(predicted,lats_for_conf,0);

%% Laterality model
lat_info = laterality_loo(hemi_lr,hemi_lr_spikes,soz_lats);
end

%% Get average L-R connectivity for each patient
%{
lr_conn = nan(npts,1);
lr_spikes = nan(npts,2);
for ip = 1:npts
    assert(abs(nanmean(mt_atlas(left,right,ip),'all') - nanmean(mt_atlas(right,left,ip),'all')) < 1e-3 ...
        || isnan(nanmean(mt_atlas(left,right,ip),'all')))
    lr_conn(ip) = nanmean(mt_atlas(left,right,ip),'all');
    lr_spikes(ip,:) = [nanmean(mt_spikes(left,ip)) nanmean(mt_spikes(right,ip))];
end
%}

if 0
    figure % odd, appears that connectivity between L MT region and R MT region higher in L sided epilepsy than in R
    plot(0.05*randn(sum(strcmp(soz_lats,'left')),1)+1,lr_conn(strcmp(soz_lats,'left')),'o')
    hold on
    plot(0.05*randn(sum(strcmp(soz_lats,'right')),1)+2,lr_conn(strcmp(soz_lats,'right')),'o')
    plot(0.05*randn(sum(strcmp(soz_lats,'bilateral')),1)+3,lr_conn(strcmp(soz_lats,'bilateral')),'o')
    
end


%% Output data for plots

nout.soz_lats = soz_lats;
nout.all_bilateral = all_bilateral;
nout.left = left;
nout.right = right;
nout.neither_lat = neither_lat;
nout.soz_non_soz_ordered_atlas = soz_non_soz_ordered_atlas;
nout.soz_all = soz_all;
nout.hemi = hemi;
nout.soz_intra = soz_intra;
nout.soz_coh_all = soz_coh_all;
nout.which_atlas = which_atlas;
nout.conf_out_fc = conf_out_fc;
nout.conf_out_spikes = conf_out_spikes;
nout.bin_soz = bin_soz;
nout.all_bin_contra_soz = all_bin_contra_soz;
nout.lat_info = lat_info;
nout.atlas_names = names;

if force_same_num == 0
if no_symm_cov
    save([plot_folder,'no_symm_cov_',which_atlas,'.mat'],'nout');
else
    save([plot_folder,'symm_cov_',which_atlas,'.mat'],'nout');
end
end

if do_plots
    
%% Compare lr conn between unilat and bilateral
if 0
    figure
    plot(1+randn(sum(unilat),1)*0.05,lr_conn(unilat),'o','linewidth',2)
    hold on
    plot(2+randn(sum(bilat),1)*0.05,lr_conn(bilat),'o','linewidth',2)
    xticks([1 2])
    xticklabels({'Unilateral','bilateral'})
    ylabel('Left mesial temporal - right mesial temporal connectivity')
    p = ranksum(lr_conn(bilat),lr_conn(unilat));
    title(sprintf('%s',get_p_text(p)))
    set(gca,'fontsize',15)

end


%% Number with left laterality and right laterality
fprintf(['\n%d patients (%1.1f%%) had left sided seizure onset, %d (%1.1f%%) '...
    'had right sided onset, %d (%1.1f%%) had bilateral onsets, and %d (%1.1f%%) '...
    'had no seizures.\n'],sum(strcmp(soz_lats,'left')),sum(strcmp(soz_lats,'left'))/length(soz_lats)*100,...
    sum(strcmp(soz_lats,'right')),sum(strcmp(soz_lats,'right'))/length(soz_lats)*100,...
    sum(strcmp(soz_lats,'bilateral')|strcmp(soz_lats,'diffuse')),sum(strcmp(soz_lats,'bilateral')|strcmp(soz_lats,'diffuse'))/length(soz_lats)*100,...
    sum(cellfun(@isempty,soz_lats)),sum(cellfun(@isempty,soz_lats))/length(soz_lats)*100);


%% Main plot
figure
set(gcf,'position',[-2023         710        1781         900])
tiledlayout(2,6,'tilespacing','compact','padding','tight')

nexttile([1 3])
pretty_matrix(all_bilateral(~neither_lat,:),...
    {'SOZ','non-SOZ'},sum(left),[],1);
colormap(gca,[0.5,0.5,0.5;1 1 1])
title('Regions with symmetric coverage')
xlabel('Patient')
ylabel('Region')

nexttile([1 3])
pretty_matrix(nanmean(soz_non_soz_ordered_atlas(~neither_lat,~neither_lat,:),3),...
    {'SOZ','non-SOZ'},sum(left),'r',0)
title('Average connectivity (symmetric coverage only)')

nexttile([1 2])
plot_type = 'scatter';
plot_paired_data(soz_all',{'SOZ','contralateral region','contralateral region'},'Average connectivity','paired',plot_type);
title('Average connectivity to SOZ vs contralateral region')

nexttile([1 2])
plot_paired_data(hemi',{'SOZ side','non-SOZ side','non-SOZ side'},'Intra-hemispheric connectivity','paired',plot_type);
title('Intra-hemispheric connectivity on side of SOZ vs non-SOZ')
nexttile([1 2])
plot_paired_data(soz_intra',{'SOZ','contralateral region','contralateral region'},'Intrinsic connectivity','paired',plot_type);
title('Intrinsic connectivity in SOZ vs contralateral region')


print(gcf,[plot_folder,'symm_',which_atlas],'-dpng')

%% coherence
figure
set(gcf,'position',[1 100 1450 350])
t=tiledlayout(1,5,'tilespacing','compact','padding','tight');

for f = 1:nfreqs
    nexttile
    paired_plot((squeeze(soz_coh_all(:,f,:))),'Coherence',{'SOZ','Contralateral region','SOZ'});
    title(freqs{f})
end
title(t,'Average coherence to SOZ vs contralateral region','fontsize',18,'fontweight','bold')
print(gcf,[plot_folder,'symm_coh_',which_atlas],'-dpng')




%% LR to check
figure
set(gcf,'position',[1 100 1140 400])
t = tiledlayout(1,3,'tilespacing','compact','padding','tight');
nexttile
plot_paired_data(hemi_lr',{'left','right','right'},'Intra-hemispheric connectivity','paired',plot_type);
title('Intra-hemispheric connectivity')

nexttile
plot_paired_data(lr_soz_all',{'left','right','right'},'Average connectivity','paired',plot_type);
title('Average connectivity to SOZ localization')

nexttile
plot_paired_data(lr_soz_intra',{'left','right','right'},'Intrinsic connectivity','paired',plot_type);
title('Intrinsic connectivity in SOZ localization')
title(t,'Connectivity left vs right','fontsize',18,'fontweight','bold')

print(gcf,[plot_folder,'symm_lr_',which_atlas],'-dpng')

%% COnfusion matrix to lateralize epilepsy
figure
set(gcf,'position',[1 100 800 370])
t = tiledlayout(1,2,'tilespacing','compact','padding','tight');

for i = 1:2 % connectivity then spikes
    
    if i == 1
        conf_out = conf_out_fc;
        ttext = 'Connectivity';
    else
        conf_out = conf_out_spikes;
        ttext = 'Spikes';
    end
    
    nexttile
    turn_nans_gray([1 0;0 1])
    colormap(gca,[0.8500, 0.3250, 0.0980;0, 0.4470, 0.7410])
    xticks(1:conf_out.nclasses)
    xticklabels(conf_out.classes)
    yticks(1:conf_out.nclasses)
    yticklabels(conf_out.classes)
    xlabel(conf_out.xlabel)
    ylabel(conf_out.ylabel)
    hold on
    for ic = 1:conf_out.nclasses
        for jc = 1:conf_out.nclasses
            text(ic,jc,sprintf('%d',conf_out.mat(jc,ic)),'horizontalalignment','center','fontsize',25,'fontweight','bold')
        end
    end
    title(sprintf('%s\nAccuracy: %1.1f%%, PPV: %1.1f%%, NPV: %1.1f%%',ttext,...
        conf_out.accuracy*100,...
        conf_out.ppv*100,conf_out.npv*100))
    set(gca,'fontsize',15)
end

print(gcf,[locations.paper_plot_folder,'symm_pred_',which_atlas],'-dpng')


end


end