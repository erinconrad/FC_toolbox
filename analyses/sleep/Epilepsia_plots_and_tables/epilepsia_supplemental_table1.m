function epilepsia_supplemental_table1

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/epilepsia/'];
int_folder = [results_folder,'analysis/backup_intermediate_Feb26_good_spikes/'];
%int_folder = [results_folder,'analysis/intermediate/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
out_folder1 = [scripts_folder,'analyses/sleep/data/'];

%% Load out file
out = load([out_folder1,'out.mat']);
out = out.out;
elec_locs = out.circ_out.all_elec_locs;
elec_lats = out.circ_out.all_elec_lats;
summ = count_locs_and_lats(elec_locs,elec_lats);

%% Get numbers
% Total number of patients with localizations
npts = sum(summ.empty_locs_lats(:,1)==0);
num_unilateral = sum(summ.loc_lat_count(:,2)==1);
num_bilateral = sum(summ.loc_lat_count(:,2)==2);
num_regions = summ.loc_lat_count(:,1);
num_perc_left = summ.lat_nums(1,:);
num_perc_right = summ.lat_nums(2,:);
num_perc_neo = summ.loc_nums(1,:);
num_perc_mt = summ.loc_nums(2,:);
num_perc_other = summ.loc_nums(3,:);

%% Prep strings
total_str = {'Number of patients: N',sprintf('%d',npts)};
lat_str = {'Implant laterality',''};
num_unilat_str = {'Unilateral implant: N (%)',sprintf('%d (%1.1f%%)',num_unilateral,...
    num_unilateral/npts*100)};
num_bilat_str = {'Bilateral implant: N (%)',sprintf('%d (%1.1f%%)',num_bilateral,...
    num_bilateral/npts*100)};
num_left_str = {'Left-sided coverage: N (%)',sprintf('%d (%1.1f%%)',...
    num_perc_left(1),num_perc_left(2)*100)};
num_right_str = {'Right-sided coverage: N (%)',sprintf('%d (%1.1f%%)',...
    num_perc_right(1),num_perc_right(2)*100)};
reg_str = {'Regional coverage',''};
num_regions_str = {'Number of regions implanted: mean (range)',sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmean(num_regions),min(num_regions),max(num_regions))};
num_neo_str = {'Temporal neocortical coverage: N (%)',sprintf('%d (%1.1f%%)',...
    num_perc_neo(1),num_perc_neo(2)*100)};
num_mt_str = {'Mesial temporal coverage: N (%)',sprintf('%d (%1.1f%%)',...
    num_perc_mt(1),num_perc_mt(2)*100)};
num_other_str = {'Other cortical coverage: N (%)',sprintf('%d (%1.1f%%)',...
    num_perc_other(1),num_perc_other(2)*100)};


all = [total_str;...
    lat_str;...
    num_unilat_str;...
    num_bilat_str;...
    num_left_str;...
    num_right_str;...
    reg_str;...
    num_regions_str;...
    num_neo_str;...
    num_mt_str;...
    num_other_str];

T = cell2table(all);
writetable(T,[out_folder,'Supplemental Table 1.csv']);

end