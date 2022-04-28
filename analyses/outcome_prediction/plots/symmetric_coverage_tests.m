function nout = symmetric_coverage_tests(which_atlas)

%{
This performs the analysis comparing connectivity to the SOZ to the region
contralateral to the SOZ, first restricting regions to only those with
bilateral coverage.
%}

%% Parameters
do_plots = 0;
freqs = {'Delta (0.5-4 Hz)','Theta (4-8 Hz)','Alpha (8-12 Hz)','Beta (12-30 Hz)','Gamma (30-80 Hz)'};

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];

bct_folder= locations.bct;
out_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];
atlas_folder = [results_folder,'analysis/atlas/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end
if ~exist(plot_folder,'dir'), mkdir(plot_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Load out file with functional connectivity and spikes as well as SOZ info
out = load([out_folder,'main_out.mat']);
out = out.out;

%% Get stuff
nfreqs = length(freqs);
rate = out.all_spikes;
soz = out.all_soz_bin;
npts = length(soz);
labels = out.all_labels;
fc = out.all_fc;
coh = out.all_coh;

% Soz lats
soz_lats = out.all_soz_lats;
right_lat = strcmp(soz_lats,'right');
left_lat = strcmp(soz_lats,'left');
%bilat = cellfun(@(x) contains(x,'bilateral') || contains(x,'diffuse'),soz_lats);
%unilat = cellfun(@(x) contains(x,'left') || contains(x,'right'),soz_lats);


%% Load atlas file
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;

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
[atlas,spikes,bin_soz,coh] = rebuild_atlas(fc,rate,atlas_elec_labels,...
    atlas_elec_regions,atlas_nums,labels,soz,coh);

%% Re-order atlas to be left then right then neither
lr_order = reorder_lr(locs,lats);

left = left(lr_order);
right = right(lr_order);
neither_lat = neither_lat(lr_order);
atlas = atlas(lr_order,lr_order,:);
coh = coh(lr_order,lr_order,:,:);
names = names(lr_order);
locs = locs(lr_order);
lats = lats(lr_order);
spikes = spikes(lr_order,:);
bin_soz = bin_soz(lr_order,:);
atlas_nums = atlas_nums(lr_order);

%% index of contralateral region for each region
contra_index = nan(size(left));
contra_index(1:sum(left)) = ([1:sum(left)])'+sum(left);
contra_index(sum(left)+1:sum(left)*2) = ([1:sum(left)])';

% make sure first half locs are same as last half locs
assert(isequal(locs(1:sum(left)),locs(sum(left)+1:sum(left)*2)))
assert(isequal(locs(contra_index(1:sum(left))),locs(contra_index(sum(left)+1:sum(left)*2))))

%% Define names corresponding to mesial temporal
switch which_atlas
    case 'aal_bernabei'
        mt_names = {'Hippocampus','Amygdala'};
    case 'brainnetome'
        mt_names = {'Amyg','Hipp'};
end

mt = contains(names,mt_names);

%% First, build symmetric coverage atlas
% Build atlas for correlation and get indices of bilateral coverage regions
[symm_cov_atlas,all_bilateral] = build_symmetric_coverage_atlas(atlas,locs,lats);
atlas = symm_cov_atlas;

% Build symmetric coverage coherence matrix
symm_cov_coh = nan(size(coh));
for ip = 1:npts
    curr_bilateral = logical(all_bilateral(:,ip));
    curr_coh = coh(:,:,:,ip);
    curr_coh(~curr_bilateral,:,:) = nan;
    curr_coh(:,~curr_bilateral,:) = nan;
    symm_cov_coh(:,:,:,ip) = curr_coh;
end
coh = symm_cov_coh;

% Build symmetric coverage spike matrix
symm_spikes = nan(size(spikes));
for ip = 1:npts
    curr_bilateral = logical(all_bilateral(:,ip));
    curr_spikes = spikes(:,ip);
    curr_spikes(~curr_bilateral) = nan;
    symm_spikes(:,ip) = curr_spikes;
end
spikes = symm_spikes;

%% Double check symmetric coverage
% Looks good!
for ip = 1:npts
    curr = atlas(:,:,ip);
    avg_rows = nanmean(curr,2);
    
    % find non nan
    non_nan = find(~isnan(avg_rows));
    
    % confirm contralateral is non nan
    for i = 1:length(non_nan)
        assert(~isnan(avg_rows(contra_index(non_nan(i))))|length(non_nan)==3) % latter is an edge case
    end
    
end

%% Build atlas ONLY containing mesial temporal regions
mt_atlas = atlas;
mt_atlas(~mt,:,:) = nan; mt_atlas(:,~mt,:) = nan;
mt_spikes = spikes;
mt_spikes(~mt,:) = nan;
mt_names = names;
mt_names(~mt) = {''};

%% Get average connectivity of SOZ (to other things) and that of contralateral region
% INitialize things
soz_intra = nan(npts,2);
lr_soz_intra = nan(npts,2);
soz_all = nan(npts,2);
lr_soz_all = nan(npts,2);
hemi = nan(npts,2);
hemi_lr = nan(npts,2);
soz_coh_all = nan(npts,nfreqs,2);
all_bin_contra_soz = zeros(size(bin_soz));

for ip = 1:npts
    
    %% SOZ connectivity
    % get soz indices
    curr_soz = bin_soz(:,ip);
    
    % get the regions contalateral to soz
    contra_soz = contra_index(curr_soz);
    contra_soz(isnan(contra_soz)) = [];
    bin_contra_soz = zeros(length(curr_soz),1);
    bin_contra_soz(contra_soz) = 1;
    bin_contra_soz = logical(bin_contra_soz);
    all_bin_contra_soz(:,ip) = bin_contra_soz;
    
    % confirm that difference in number of soz and contra soz is just those
    % regions that have no contralateral thing 
    switch which_atlas
        case 'aal_bernabei'
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
end

%% Build a SOZ - non SOZ laterality ordered atlas
% This is for plotting purposes
soz_non_soz_ordered_atlas = build_soz_ordered_atlas(atlas,left,right,right_lat,left_lat);

%% get confusion matrix for connectivity
hemi_diff = hemi_lr(:,1) - hemi_lr(:,2);
predicted = cell(length(soz_lats),1);
predicted(hemi_diff > 0) = {'right'};
predicted(hemi_diff < 0) = {'left'};
empty = hemi_diff == 0 | isnan(hemi_diff) | (~strcmp(soz_lats,'right') & ~strcmp(soz_lats,'left'));
predicted(empty) = [];
lats_for_conf = soz_lats;
lats_for_conf(empty) = [];

if 0
    table(lats_for_conf,predicted)
end

conf_out_fc = confusion_matrix(predicted,lats_for_conf,0);

%% get confusion matrix for spikes
% Note that this is also using a symmetric coverage map!
hemi_lr_spikes = [nanmean(spikes(left,:),1)',nanmean(spikes(right,:),1)'];
hemi_diff = hemi_lr_spikes(:,1) - hemi_lr_spikes(:,2);
predicted = cell(length(soz_lats),1);
predicted(hemi_diff < 0) = {'right'};
predicted(hemi_diff > 0) = {'left'};
empty = hemi_diff == 0 | isnan(hemi_diff) | (~strcmp(soz_lats,'right') & ~strcmp(soz_lats,'left'));
predicted(empty) = [];
lats_for_conf = soz_lats;
lats_for_conf(empty) = [];

if 0
    table(lats_for_conf,predicted)
end

conf_out_spikes = confusion_matrix(predicted,lats_for_conf,0);

%% Get average L-R connectivity for each patient
lr_conn = nan(npts,1);
lr_spikes = nan(npts,2);
for ip = 1:npts
    assert(abs(nanmean(mt_atlas(left,right,ip),'all') - nanmean(mt_atlas(right,left,ip),'all')) < 1e-3 ...
        || isnan(nanmean(mt_atlas(left,right,ip),'all')))
    lr_conn(ip) = nanmean(mt_atlas(left,right,ip),'all');
    lr_spikes(ip,:) = [nanmean(mt_spikes(left,ip)) nanmean(mt_spikes(right,ip))];
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
save([plot_folder,'symm_cov_',which_atlas,'.mat'],'nout');

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

print(gcf,[plot_folder,'symm_pred_',which_atlas],'-dpng')

%% Output text
fprintf(fid,'<p><b>Changes in spikes with seizures</b>');

end


end