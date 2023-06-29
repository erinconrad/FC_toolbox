function make_table_1

%% Parameters
rm_non_temporal = 1;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Initialize results file
fname = [plot_folder,'results.html'];
fid = fopen(fname,'a');

%% Run the lr_mt to extract features
[T,features,Ts] =  lr_mt(3,0);
empty_class = cellfun(@isempty,T.soz_lats);
T(empty_class,:) = [];

%% Load the pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Go through and remove non-temporal patients
soz_loc = T.soz_locs;
tle = contains(soz_loc,'temporal');
etle = ~tle;

if rm_non_temporal
    % Remove from T (all subsequent things map to the patients in T)
    T(etle,:) = [];
    tle(etle,:) = [];
    etle(etle,:) = [];

end

%% Grab demographic variables from table
% Get the patient names
name = T.names;
npts = length(name);

% HUP vs musc
site = cell(npts,1);
is_hup = contains(name,'HUP');
is_musc = contains(name,'MP');
site(is_hup) = {'HUP'};
site(is_musc) = {'MUSC'};
nhup = sum(is_hup);
nmusc = sum(is_musc);

% Get outcomes
engel_yr1 = cellfun(@(x) parse_outcome_num(x,'engel'),T.engel_yr1);
engel_yr2 = cellfun(@(x) parse_outcome_num(x,'engel'),T.engel_yr2);
ilae_yr1 = cellfun(@(x) parse_outcome_num(x,'ilae'),T.ilae_yr1);
ilae_yr2 = cellfun(@(x) parse_outcome_num(x,'ilae'),T.ilae_yr2);

% Get surgery and soz locs and lats
surg = T.surgery;
resection = contains(surg,'Resection');
ablation = contains(surg,'ablation');
device = contains(surg,'RNS') | contains(surg,'DBS') | contains(surg,'VNS');


soz_lat = T.soz_lats;

left = strcmp(soz_lat,'left');
right = strcmp(soz_lat,'right');
bilateral = strcmp(soz_lat,'bilateral');
assert(sum(left)+sum(right)+sum(bilateral)==npts)

% nummber of symmetric mt electrodes
n_symmetric = T.n_symmetric;

% ADD IN PERCENT OF TIMES CLASSIFIED AS ASLEEP
n_wake = T.n_wake;
n_sleep = T.n_sleep;
n_connected = T.n_connected;
perc_wake = n_wake./n_connected;
perc_sleep = n_sleep./n_connected;


%% Get the other demographic variables from the pt.mat
% prep them
female = nan(npts,1);
age_onset = nan(npts,1);
age_implant = nan(npts,1);
has_dkt = zeros(npts,1);
has_atropos = zeros(npts,1);
rid = nan(npts,1);


% loop over patients
for ip = 1:length(pt)
    % see if the name matches any of the main names
    ip_name = pt(ip).name;
    r = strcmp(ip_name,name);

    if sum(r) ~= 1, continue; end

    if isfield(pt(ip),'atropos')
        if ~isempty(pt(ip).atropos)
            has_atropos(r) = 1;
        end
    end

    if isfield(pt(ip),'dkt')
        if ~isempty(pt(ip).dkt)
            has_dkt(r) = 1;
        end
    end

    if ~isempty(pt(ip).rid)
        rid(r) = pt(ip).rid;
    end

    % get the demographics
    if ~isfield(pt(ip),'clinical') || isempty(pt(ip).clinical)
        continue;
    end

    sex = pt(ip).clinical.sex;
    if strcmp(sex,'Female') == 1
        female(r) = 1;
    elseif strcmp(sex,'Male') == 1
        female(r) = 0;
    else
        female(r) = nan;
    end

    age_onset(r) = pt(ip).clinical.age_onset;
    age_implant(r) = pt(ip).clinical.age_implant;

    

end

%% Maybe get some preimplant data from the manual validation file
% MRI, scalp seizure laterality, scalp spike laterality, PET
mT = readtable('Manual validation.xlsx','Sheet','Pre-implant data');
no1_lat = cell(npts,1);
no2_lat = cell(npts,1);
for i = 1:npts
    % See if you can find the patient in this table
    curr_name = name{i};

    r = strcmp(curr_name,mT.name);

    if r ~=1, continue; end

    % Get the laterality hypotheses
    no1_lat{i} = lower(mT.x_1PreimplantHypothesisLaterality_left_Right_Bilateral_Broad_NA{r});
    no2_lat{i} = lower(mT.x_2PreimplantHypothesisLaterality_left_Right_Bilateral_Broad_NA{r});

    

end
no1_lat(cellfun(@isempty,no1_lat)) = {''};
no2_lat(cellfun(@isempty,no2_lat)) = {''};

% Do some logic to decide if the patient's top two hypotheses contain
% either 1) a bilateral hypothesis or 2) discordant lateralities
bilat_hypothesis = ismember(no1_lat,{'bilateral','broad'}) | ismember(no2_lat,{'bilateral','broad'});
bilat_hypothesis = double(bilat_hypothesis);
bilat_hypothesis(strcmp(no1_lat,'') & strcmp(no2_lat,'')) = nan;
discordant_hypotheses = cellfun(@(x,y) ~strcmp(x,y),no1_lat,no2_lat);
bilat_or_discordant = bilat_hypothesis == 1 | discordant_hypotheses == 1;
bilat_or_discordant = double(bilat_or_discordant);
bilat_or_discordant(strcmp(no1_lat,'') & strcmp(no2_lat,'')) = nan;

if 0
    table(name,no1_lat,no2_lat,bilat_or_discordant)
end

%% Put the table together
% Planning to have 3 total columns: the first column says the thing, the
% second is the data for HUP, the third is the data for MUSC
total_str = {'Total: N',sprintf('%d',nhup),sprintf('%d',nmusc)};
female_str = {'Female: N (%)',sprintf('%d (%1.1f%%)',sum(female==1 & is_hup),sum(female==1 & is_hup)/sum(~isnan(female(is_hup)))*100),...
    sprintf('%d (%1.1f%%)',sum(female==1 & is_musc),sum(female==1 & is_musc)/sum(~isnan(female(is_musc)))*100)};
age_onset_str = {'Age at onset in years: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(age_onset(is_hup)),min(age_onset(is_hup)),max(age_onset(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(age_onset(is_musc)),min(age_onset(is_musc)),max(age_onset(is_musc)))};
age_implant_str = {'Age at implant in years: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(age_implant(is_hup)),min(age_implant(is_hup)),max(age_implant(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(age_implant(is_musc)),min(age_implant(is_musc)),max(age_implant(is_musc)))};
n_discordant_str = {'Bilateral or discordant pre-implant hypotheses: N (%)',...
    sprintf('%d (%1.1f%%)',sum(bilat_or_discordant==1 & is_hup),sum(bilat_or_discordant==1 & is_hup)/sum(~isnan(bilat_or_discordant(is_hup)))*100),...
    sprintf('%d (%1.1f%%)',sum(bilat_or_discordant==1 & is_musc),sum(bilat_or_discordant==1 & is_musc)/sum(~isnan(bilat_or_discordant(is_musc)))*100)};
n_elecs_str = {'Symmetric mesial temporal-targeted contacts: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(n_symmetric(is_hup)),min(n_symmetric(is_hup)),max(n_symmetric(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(n_symmetric(is_musc)),min(n_symmetric(is_musc)),max(n_symmetric(is_musc)))};
loc_str = {'SOZ localization (clinician determination)','',''};
tle_str = {'Temporal: N (%)',...
    sprintf('%d (%1.1f%%)',sum(tle==1 & is_hup),sum(tle==1 & is_hup)/sum(is_hup)*100),...
    sprintf('%d (%1.1f%%)',sum(tle==1 & is_musc),sum(tle==1 & is_musc)/sum(is_musc)*100)};
etle_str = {'Extratemporal: N (%)',...
    sprintf('%d (%1.1f%%)',sum(etle==1 & is_hup),sum(etle==1 & is_hup)/sum((is_hup))*100),...
    sprintf('%d (%1.1f%%)',sum(etle==1 & is_musc),sum(etle==1 & is_musc)/sum(is_musc)*100)};
lat_str = {'SOZ lateralization (clinician determination)','',''};
left_str = {'Left: N (%)',...
    sprintf('%d (%1.1f%%)',sum(left==1 & is_hup),sum(left==1 & is_hup)/sum(~isnan(left(is_hup)))*100),...
    sprintf('%d (%1.1f%%)',sum(left==1 & is_musc),sum(left==1 & is_musc)/sum(~isnan(left(is_musc)))*100)};
right_str = {'Right: N (%)',...
    sprintf('%d (%1.1f%%)',sum(right==1 & is_hup),sum(right==1 & is_hup)/sum(~isnan(right(is_hup)))*100),...
    sprintf('%d (%1.1f%%)',sum(right==1 & is_musc),sum(right==1 & is_musc)/sum(~isnan(right(is_musc)))*100)};
bilat_str = {'Bilateral: N (%)',...
    sprintf('%d (%1.1f%%)',sum(bilateral==1 & is_hup),sum(bilateral==1 & is_hup)/sum((is_hup))*100),...
    sprintf('%d (%1.1f%%)',sum(bilateral==1 & is_musc),sum(bilateral==1 & is_musc)/sum((is_musc))*100)};
surg_str = {'Surgery performed','',''};
resection_str = {'Resection: N (%)',...
    sprintf('%d (%1.1f%%)',sum(resection==1 & is_hup),sum(resection==1 & is_hup)/sum(is_hup)*100),...
    sprintf('%d (%1.1f%%)',sum(resection==1 & is_musc),sum(resection==1 & is_musc)/sum(is_musc)*100)};
ablation_str = {'Ablation: N (%)',...
    sprintf('%d (%1.1f%%)',sum(ablation==1 & is_hup),sum(ablation==1 & is_hup)/sum(is_hup)*100),...
    sprintf('%d (%1.1f%%)',sum(ablation==1 & is_musc),sum(ablation==1 & is_musc)/sum(is_musc)*100)};
device_str = {'Device: N (%)',...
    sprintf('%d (%1.1f%%)',sum(device==1 & is_hup),sum(device==1 & is_hup)/sum(is_hup)*100),...
    sprintf('%d (%1.1f%%)',sum(device==1 & is_musc),sum(device==1 & is_musc)/sum(is_musc)*100)};
engel_str = {'Engel outcome','',''};
engel_one_str = {'Year 1: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(engel_yr1(is_hup)),min(engel_yr1(is_hup)),max(engel_yr1(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(engel_yr1(is_musc)),min(engel_yr1(is_musc)),max(engel_yr1(is_musc)))};
engel_two_str = {'Year 2: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(engel_yr2(is_hup)),min(engel_yr2(is_hup)),max(engel_yr2(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(engel_yr2(is_musc)),min(engel_yr2(is_musc)),max(engel_yr2(is_musc)))};
ilae_str = {'ILAE outcome','',''};
ilae_one_str = {'Year 1: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(ilae_yr1(is_hup)),min(ilae_yr1(is_hup)),max(ilae_yr1(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(ilae_yr1(is_musc)),min(ilae_yr1(is_musc)),max(ilae_yr1(is_musc)))};
ilae_two_str = {'Year 2: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(ilae_yr2(is_hup)),min(ilae_yr2(is_hup)),max(ilae_yr2(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(ilae_yr2(is_musc)),min(ilae_yr2(is_musc)),max(ilae_yr2(is_musc)))};

all = [total_str;...
    female_str;...
    age_onset_str;...
    age_implant_str;...
    n_discordant_str;...
    n_elecs_str;...
    loc_str;...
    tle_str;...
    etle_str;...
    lat_str;...
    left_str;...
    right_str;...
    bilat_str;...
    surg_str;...
    resection_str;...
    ablation_str;...
    device_str;...
    engel_str;...
    engel_one_str;...
    engel_two_str;...
    ilae_str;...
    ilae_one_str;...
    ilae_two_str];

T2 = cell2table(all);
writetable(T2,[plot_folder,'Table1.csv']);

%% Get inclusion/exclusion numbers
% How many total HUP patients did I look at?
n_hup_total = sum(contains(Ts.names,'HUP'));
hup_pts = contains(Ts.names,'HUP');

% musc
musc_pts = contains(Ts.names,'MP');
n_musc_total = sum(contains(Ts.names,'MP'));

% How many did I exclude due to no symmetric mesial temporal contacts?
n_exclude_no_symm_hup = sum(hup_pts & Ts.n_symmetric == 0);
n_exclude_no_symm_musc = sum(musc_pts & Ts.n_symmetric == 0);

% How many did I exclude due to etle
n_etle_hup = sum(hup_pts & ~contains(Ts.soz_locs,'temporal') & (Ts.n_symmetric ~= 0));
n_etle_musc = sum(musc_pts & ~contains(Ts.soz_locs,'temporal') & (Ts.n_symmetric ~= 0));

% How many did I exclude due to most of the EEG disconnected?
n_exclude_no_sleep_conn_hup = sum(hup_pts & Ts.most_disconnected); % 0
n_exclude_no_sleep_conn_musc = sum(musc_pts & Ts.most_disconnected & Ts.n_symmetric ~= 0 & contains(Ts.soz_locs,'temporal'));
assert(n_exclude_no_sleep_conn_hup == 0)
assert(n_exclude_no_sleep_conn_musc == 1)

% total I think I should have excluded
n_hup_excluded = n_exclude_no_symm_hup + n_etle_hup;
n_hup_remaining = n_hup_total - n_hup_excluded;
n_musc_excluded = n_exclude_no_symm_musc + n_etle_musc + n_exclude_no_sleep_conn_musc;
n_musc_remaining = n_musc_total - n_musc_excluded;

% Make sure it adds up to expected number of included patients
assert(sum(contains(T.names,'HUP'))==n_hup_remaining)
assert(sum(contains(T.names,'MP'))==n_musc_remaining)

%% Summary results
fprintf(fid,'<p><br><b>Clinical information and summary of intracranial recording</b></br>');
fprintf(fid,['We examined %d consecutive patients at HUP for our univariate analyses and model '...
    'development and internal validation. Of these, we excluded %d patients who did not have symmetric bitemporal '...
    'electrode coverage, and an additional %d patients who were clinically determined to have '...
    'extratemporal seizure onset zones, yielding %d HUP patients included for analysis.'],...
    n_hup_total,n_exclude_no_symm_hup,n_etle_hup,n_hup_remaining);

fprintf(fid,[' We examined %d consecutive patients at MUSC for external model validation. We excluded %d patients '...
    'with no symmetric bitemporal electrode coverage, an additional %d patient with extra-temporal seizure onset, '...
    'and an additional %d patient for whom most of the EEG was disconnected, yielding %d remaining MUSC patients who '...
    'met inclusion criteria (Table 1).'],...
    n_musc_total,n_exclude_no_symm_musc,n_etle_musc,n_exclude_no_sleep_conn_musc,n_musc_remaining);

fprintf(fid,[' As described in Table 1, the majority of patients underwent bilateral electrode implantation '...
    ' either because one of the two primary pre-implantation hypotheses was bilateral, '...
    'or because the two primary pre-implant hypotheses were discordant lateralities (e.g., '...
    'the first hypothesis was left-sided and the second hypothesis was right-sided).</p>']);

%% COmpare duration between left and right
if 0
    duration = age_implant-age_onset;
    unpaired_plot(duration(left),duration(right),{'left','right'},'duration')
    % looks like no difference
end

%% Make another table with missing atlas data
mT = table(rid,name,has_atropos,has_dkt);
writetable(mT,[plot_folder,'missingAtlasTable.csv']);

%% Prelim data for R01 aims
a = cellfun(@(x) regexp(x,'\d*','Match'),name);
b = cellfun(@(x) str2num(x), a);
N = sum(contains(name,'HUP')& ((b>=159 & b<=199) | (b == 140 | b ==143))); % 2018-2019
Nbilat = sum((contains(name,'HUP')& ((b>=159 & b<=199) | (b == 140 | b ==143))) & bilateral == 1);
%{
Looks like we started doing mostly stereo with HUP127, November 2016. Going
5 years out, this would take us up to Nov 2021, or HUP226. There are 57
between hUP127 and HUP224, which is as far as I processed. HUP225 did not
have bilateral MT, but HUP226 did. So 58 over 5 years. 11.6/year. 23 in 2
years. However, 37% of these were bilateral. So only 14.5 unilateral in 2
years.......

Ok what if I only look at 2018 and 2019 (two normal years, assuming 2020,
and 2021 screwed up due to COVID). HUP159-199 + 2 stragglers (HUP140 and
143). Then I get 29 bilateral MT  implants, 9 of whom had bilateral SOZ 
(31%). So then I can estimate 20 unilateral patients with bilateral MT implants,
 which gets 90% power.
%}


end