function sw_adr_atlas(which_atlas)

%% Get file locs
locations = fc_toolbox_locs;

bct_folder= locations.bct;
atlas_folder = locations.paper_data_folder;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Load out file with functional connectivity and spikes as well as SOZ info
out = load([locations.paper_data_folder,'main_out.mat']);
out = out.out;

%% Get stuff
labels = out.all_labels;
npts = length(labels);
all_ad_wake = out.all_ad_wake;
all_ad_sleep = out.all_ad_sleep;
all_soz_bin = out.all_soz_bin;

%% Load atlas file
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;

assert(isequal(out.all_names,atlas_out.pt_names))

%% get atlas stuff (indictates which electrodes are in which atlas regions)
atlas_elec_labels = atlas_out.elecs_labels;
atlas_elec_regions = atlas_out.elecs_atlas;
atlas_nums = atlas_out.atlas_nums;
names = atlas_out.atlas_names;

%% Put AD and SOZ into atlas
[ad_wake_atlas,~] = adr_atlas(all_ad_wake,atlas_elec_labels,atlas_elec_regions,atlas_nums,labels,all_soz_bin);
[ad_sleep_atlas,soz_atlas] = adr_atlas(all_ad_sleep,atlas_elec_labels,atlas_elec_regions,atlas_nums,labels,all_soz_bin);


%% Normalize across patients
ad_wake_norm = ad_wake_atlas;
ad_sleep_norm = ad_sleep_atlas;
%ad_wake_norm = (ad_wake_atlas-nanmean(ad_wake_atlas,2))./nanstd(ad_wake_atlas,[],2);
%ad_sleep_norm = (ad_sleep_atlas-nanmean(ad_sleep_atlas,2))./nanstd(ad_sleep_atlas,[],2);

soz_non_wake = nan(npts,2);
soz_non_sleep = nan(npts,2);
for ip = 1:npts
    curr_soz = soz_atlas(:,ip);
    curr_ad_wake_norm = ad_wake_norm(:,ip);
    soz_non_wake(ip,1) = nanmean(curr_ad_wake_norm(curr_soz==1));
    soz_non_wake(ip,2) = nanmean(curr_ad_wake_norm(curr_soz==0));
    
    curr_ad_sleep_norm = ad_sleep_norm(:,ip);
    soz_non_sleep(ip,1) = nanmean(curr_ad_sleep_norm(curr_soz==1));
    soz_non_sleep(ip,2) = nanmean(curr_ad_sleep_norm(curr_soz==0));
end

figure
nexttile
paired_plot(soz_non_wake,'Normalized wake AD ratio',{'SOZ','non-SOZ'})
nexttile
paired_plot(soz_non_sleep,'Normalized sleep AD ratio',{'SOZ','non-SOZ'})

end

function [ad_atlas,soz_atlas] = adr_atlas(ad,atlas_labels,regions,atlas_nums,labels,soz)

npts = length(ad);
nregions = length(atlas_nums);

%% Remove ekg from atlas
for ip = 1:npts
    ekg = find_non_intracranial(atlas_labels{ip});
    atlas_labels{ip}(ekg) = [];
    regions{ip}(ekg) = [];
    
    assert(isequal(atlas_labels{ip},labels{ip})); % make sure electrode labels line up
end

ad_atlas = nan(nregions,npts);
soz_atlas = zeros(nregions,npts);

% loop over patients
for ip = 1:npts
    curr_ad = ad{ip};
    curr_soz = soz{ip};
    curr_regions = regions{ip};
    
    % Loop over regions
    for ir = 1:nregions
        
        % get this region num
        curr_num_i = atlas_nums(ir);
        
        % get the channels in this region
        curr_chs_i = ismember(curr_regions,curr_num_i);
        
        % soz (1 if any electrodes in that region are SOZ)
        soz_atlas(ir,ip) = any(curr_soz(curr_chs_i)==1);
        
        % spikes (mean ad rate of channels in that region)
        ad_atlas(ir,ip) = nanmean(curr_ad(curr_chs_i));
    end
    
end


end