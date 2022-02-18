function fine_localization

%{
At individual region level, get node strength. Then normalize across
patients. Do SOZ regions have higher or lower connectivity than non SOZ
regions?
%} 

%% Parameters
which_atlas = 'brainnetome';%'aal_bernabei';%'brainnetome';% %'aal';'aal_bernabei';
plot_type = 'scatter';
coverage_limit = 20;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];
bct_folder= locations.bct;
out_folder = [results_folder,'analysis/atlas/'];
if ~exist(out_folder,'dir'), mkdir(out_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Load atlas
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;

atlas = out.atlas;
names = out.atlas_names;
pt_names = out.pt_names;
atlas_nums = out.atlas_nums;
nregions = length(names);
assert(nregions==size(atlas,1))
npts = size(atlas,3);
nelecs = out.n_elecs_all;
sozs = out.sozs;
bin_soz = (cell2mat(cellfun(@(x) ismember(atlas_nums',x),sozs,'uniformoutput',false)))';


%% Indicate regions with low coverage
any_coverage = squeeze(any(~isnan(atlas),2));

% find regions with low coverage
low_coverage = sum(any_coverage,2) < coverage_limit;


if 0
    figure
    nexttile
    turn_nans_gray(any_coverage)
    
    nexttile
    turn_nans_gray(any_coverage(~low_coverage,:))
    yticks(1:sum(~low_coverage))
    yticklabels(names(~low_coverage))
end

%% Restrict atlas and soz
atlas = atlas(~low_coverage,~low_coverage,:);
bin_soz = bin_soz(~low_coverage,:);
names = names(~low_coverage);

%% Remove white matter
if 1
white_matter = ismember(names,'White_Matter');
names = names(~white_matter);
bin_soz = bin_soz(~white_matter,:);
atlas = atlas(~white_matter,~white_matter,:);
end

if 0
    figure
    nexttile
    turn_nans_gray(nanmean(atlas,3))
    
    nexttile
    turn_nans_gray(bin_soz)
end

%% Ns
ns = squeeze(nanmean(atlas,2));
z = (ns-nanmean(ns,2))./nanstd(ns,[],2);

if 0
    figure
    nexttile
    turn_nans_gray(z)
    yticks(1:size(z,1))
    yticklabels(names)
    
    nexttile
    turn_nans_gray(bin_soz)
end

%% Z soz vs not
soz_not = nan(npts,2);
soz_not_raw = nan(npts,2);
for ip = 1:npts
    curr_soz = bin_soz(:,ip);
    soz_not(ip,:) = [nanmean((z(curr_soz,ip))), nanmean((z(~curr_soz,ip)))];
    soz_not_raw(ip,:) = [nanmean((ns(curr_soz,ip))), nanmean((ns(~curr_soz,ip)))];
end

figure
nexttile
stats = plot_paired_data((soz_not)',{'SOZ','non-SOZ','non-SOZ'},'Normalized connectivity','paired',plot_type);

nexttile

stats = plot_paired_data((soz_not_raw)',{'SOZ','non-SOZ','non-SOZ'},'Raw connectivity','paired',plot_type);

end