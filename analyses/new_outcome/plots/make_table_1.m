function make_table_1

% This can't be run by the general user because it requires pt.mt

%% Parameters
rm_non_temporal = 1;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
inter_folder = locations.el_data_folder;
alfredo_folder = [locations.main_folder,'Alfredo_code/fmri_analysis_AL_3_28_23/'];
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

%% Load the file containing intermediate data
mt_data = load([inter_folder,'mt_out_epilepsy_laterality.mat']);
mt_data = mt_data.out;

%% Run the lr_mt to extract features
[T,features,Ts] =  lr_mt(mt_data,3,0);
empty_class = cellfun(@isempty,T.soz_lats);
T(empty_class,:) = [];

%% Load the pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load Alfredo file for fMRI patients
f1T = readtable([alfredo_folder,'df.csv']);

%% Load clinical file for additional info for fMRI patients
f2T = readtable([data_folder,'clinical_info/all_rids.csv']);

%% Load MUSC clinical folder
muscT = readtable([inter_folder,'LEN patient list research erin.xlsx']);

% Get musc outcomes, surg
for ip = 1:size(T,1)

    if ~contains(T.names,'MP')
        continue
    end

    % find matching musc table patient
    musc_row = contains(muscT.LENID_,T.names{ip});
    if sum(musc_row) == 1

        % get surgery
        if strcmp(muscT.TypeOfSurgery{musc_row},'ATL') || contains(muscT.TypeOfSurgery{musc_row},'Resection')
            T.surgery{ip} = 'Resection';
        elseif strcmp(muscT.TypeOfSurgery{musc_row},'LITT') || contains(muscT.TypeOfSurgery{musc_row},'ablation')
            T.surgery{ip} = 'Laser ablation';
        elseif strcmp(muscT.TypeOfSurgery{musc_row},'RNS')
            T.surgery{ip} = 'RNS';
        end

        if contains(muscT.EngelYear1{musc_row},'N/a') % don't change format if empty outcome
            T.engel_yr1{ip} = '';
            T.engel_yr2{ip} = '';
            T.ilae_yr1{ip} = '';
            T.ilae_yr2{ip} = '';
            continue
        end

        % get outcomes
        T.engel_yr1{ip} = muscT.EngelYear1{musc_row};
        T.engel_yr2{ip} = muscT.EngelYear2{musc_row};
        T.ilae_yr1{ip} = muscT.ILAEYear1{musc_row};
        T.ilae_yr2{ip} = muscT.ILAEYear2{musc_row};

        

    end


end

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

% Get surgery and soz locs and lats
T.surgery(strcmp(T.names,'HUP224')) = {'none'}; % 224 did not have surgery, just planned to get it
surg = T.surgery;
resection = contains(surg,'Resection');
ablation = contains(surg,'ablation');
device = contains(surg,'RNS') | contains(surg,'DBS') | contains(surg,'VNS');


% Get outcomes
engel_yr1 = cellfun(@(x) parse_outcome_num(x,'engel'),T.engel_yr1);
engel_yr2 = cellfun(@(x) parse_outcome_num(x,'engel'),T.engel_yr2);
ilae_yr1 = cellfun(@(x) parse_outcome_num(x,'ilae'),T.ilae_yr1);
ilae_yr2 = cellfun(@(x) parse_outcome_num(x,'ilae'),T.ilae_yr2);



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

%% Get MUSC demographics
for ip = 1:size(muscT,1)
    musc_name = muscT.LENID_{ip};
    part_name = musc_name(4:end);

    % see if the name matches any of the main names
    r = strcmp(part_name,name);
    if sum(r) ~= 1, continue; end

    sex = muscT.Sex{ip};
    if strcmp(sex,'F') == 1
        female(r) = 1;
    elseif strcmp(sex,'M') == 1
        female(r) = 0;
    else
        error('what')
    end

    age_onset(r) = muscT.AgeOfSeizureOnset(ip);
    age_implant(r) = muscT.AgeAtImplant(ip);
    
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

    if sum(r) ~=1, continue; end

    % Get the laterality hypotheses
    no1_lat{i} = lower(mT.x_1PreimplantHypothesisLaterality_left_Right_Bilateral_Broad_NA{r});
    no2_lat{i} = lower(mT.x_2PreimplantHypothesisLaterality_left_Right_Bilateral_Broad_NA{r});

    

end

%% Add preimplant hypotheses from MUSC
for i = 1:npts
    % See if you can find the patient in this table
    curr_name = name{i};

    r = contains(muscT.LENID_,['3T_',curr_name]);
    if sum(r) ~=1, continue; end

    % Get the laterality hypotheses
    hyp = muscT.Primary2HypothesesForEpilepsyLocalization{r};
    C = strsplit(hyp,'2');
    if length(C) ~=2, error('what'); end

    for ic = 1:2
        if contains(C{ic},'Right')
            if ic == 1
                no1_lat{i} = 'right';
            elseif ic == 2
                no2_lat{i} = 'right';
            end
        elseif contains(C{ic},'Left')
            if ic == 1
                no1_lat{i} = 'left';
            elseif ic == 2
                no2_lat{i} = 'left';
            end
        else
            if strcmp(curr_name,'MP0029')
                no2_lat{i} = 'right';
            else
                error('what');
            end
        end
    end


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

%% Get information for the fMRI patients
fmri_non_control = strcmp(f1T.Control,'No'); nfmri_no_control = sum(fmri_non_control);
fmri_rids = f1T.Subject(fmri_non_control);
fmri_rid_nums = cellfun(@(x) str2num(x(end-2:end)),fmri_rids);

% get age
age_fmri = f1T.Age(fmri_non_control);

% age of onset
age_onset_fmri = f1T.SzOnset(fmri_non_control);

% sex
sex_fmri = f1T.Sex(fmri_non_control);

% lat
lat_fmri = f1T.Final_Lat(fmri_non_control);

% get corresponding rows in my demographics table
[Lia,Locb]=ismember(fmri_rid_nums,f2T.record_id); % Locb has the corresponding rows
assert(sum(Lia==1)==length(Lia))
assert(isequal(f2T.record_id(Locb),fmri_rid_nums))

% Get number with resection, ablation, and device
resection_fmri = (f2T.outcome_proctype___1(Locb)==1 | f2T.outcome_proctype___2(Locb)==1);
ablation_fmri = (f2T.outcome_proctype___4(Locb)==1);
device_fmri = (f2T.outcome_proctype___3(Locb)==1 | f2T.outcome_proctype___5(Locb)==1 | f2T.outcome_proctype___6(Locb)==1);
nresection_fmri = sum(resection_fmri);
nablation_fmri = sum(ablation_fmri);
ndevice_fmri = sum(device_fmri);

% outcomes
tengel1 = f2T.demog_year(Locb & (resection_fmri | ablation_fmri));
tilae1 = f2T.demog_ilae1year(Locb& (resection_fmri | ablation_fmri));
tengel2 = f2T.demog_years2(Locb& (resection_fmri | ablation_fmri));
tilae2 = f2T.demog_ilae2years(Locb& (resection_fmri | ablation_fmri));
ilae1_fmri = get_ilae(tilae1);
engel1_fmri = get_engel(tengel1);
ilae2_fmri = get_ilae(tilae2);
engel2_fmri = get_engel(tengel2);
engel_yr1_fmri = cellfun(@(x) parse_outcome_num(x,'engel'),engel1_fmri);
engel_yr2_fmri = cellfun(@(x) parse_outcome_num(x,'engel'),engel2_fmri);
ilae_yr1_fmri = cellfun(@(x) parse_outcome_num(x,'ilae'),ilae1_fmri);
ilae_yr2_fmri = cellfun(@(x) parse_outcome_num(x,'ilae'),ilae2_fmri);

if 0
    table(fmri_rids,f2T.outcome_proctype___1(Locb)==1 | f2T.outcome_proctype___2(Locb)==1,...
        f2T.outcome_proctype___4(Locb)==1,...
        f2T.outcome_proctype___3(Locb)==1 | f2T.outcome_proctype___5(Locb)==1 | f2T.outcome_proctype___6(Locb)==1,...
        ilae1_fmri,ilae2_fmri,engel1_fmri,engel2_fmri,...
        'VariableNames',{'RID','Resection','Ablation','Device','ILAE1','ILAE2','Engel1','Engel2'})
end

%% Put the table together
% Planning to have 4 total columns: the first column says the thing, the
% second is the data for HUP, the third is the data for MUSC, the fourth
% column is for fMRI
total_str = {'Total: N',sprintf('%d',nhup),sprintf('%d',nmusc),sprintf('%d',sum(fmri_non_control))};
female_str = {'Female: N (%)',sprintf('%d (%1.1f%%)',sum(female==1 & is_hup),sum(female==1 & is_hup)/sum(~isnan(female(is_hup)))*100),...
    sprintf('%d (%1.1f%%)',sum(female==1 & is_musc),sum(female==1 & is_musc)/sum(~isnan(female(is_musc)))*100),...
    sprintf('%d (%1.1f%%)',sum(strcmp(sex_fmri,'F')),sum(strcmp(sex_fmri,'F'))/sum(fmri_non_control)*100)};
age_onset_str = {'Age at onset in years: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(age_onset(is_hup)),min(age_onset(is_hup)),max(age_onset(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(age_onset(is_musc)),min(age_onset(is_musc)),max(age_onset(is_musc))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(age_onset_fmri),min(age_onset_fmri),max(age_onset_fmri))};
age_implant_str = {'Age at implant in years: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(age_implant(is_hup)),min(age_implant(is_hup)),max(age_implant(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(age_implant(is_musc)),min(age_implant(is_musc)),max(age_implant(is_musc))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(age_fmri),min(age_fmri),max(age_fmri))};
n_discordant_str = {'Bilateral or discordant pre-implant hypotheses: N (%)',...
    sprintf('%d (%1.1f%%)',sum(bilat_or_discordant==1 & is_hup),sum(bilat_or_discordant==1 & is_hup)/sum(~isnan(bilat_or_discordant(is_hup)))*100),...
    sprintf('%d (%1.1f%%)',sum(bilat_or_discordant==1 & is_musc),sum(bilat_or_discordant==1 & is_musc)/sum(~isnan(bilat_or_discordant(is_musc)))*100),...
    'NA'};
n_elecs_str = {'Symmetric mesial temporal-targeted contacts: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(n_symmetric(is_hup)),min(n_symmetric(is_hup)),max(n_symmetric(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(n_symmetric(is_musc)),min(n_symmetric(is_musc)),max(n_symmetric(is_musc))),...
    'NA'};
%{
loc_str = {'SOZ localization (clinician determination)','','',''};

tle_str = {'Temporal: N (%)',...
    sprintf('%d (%1.1f%%)',sum(tle==1 & is_hup),sum(tle==1 & is_hup)/sum(is_hup)*100),...
    sprintf('%d (%1.1f%%)',sum(tle==1 & is_musc),sum(tle==1 & is_musc)/sum(is_musc)*100)};
etle_str = {'Extratemporal: N (%)',...
    sprintf('%d (%1.1f%%)',sum(etle==1 & is_hup),sum(etle==1 & is_hup)/sum((is_hup))*100),...
    sprintf('%d (%1.1f%%)',sum(etle==1 & is_musc),sum(etle==1 & is_musc)/sum(is_musc)*100)};
%}
lat_str = {'SOZ lateralization (clinician determination)','','',''};
left_str = {'Left: N (%)',...
    sprintf('%d (%1.1f%%)',sum(left==1 & is_hup),sum(left==1 & is_hup)/sum(~isnan(left(is_hup)))*100),...
    sprintf('%d (%1.1f%%)',sum(left==1 & is_musc),sum(left==1 & is_musc)/sum(~isnan(left(is_musc)))*100),...
    sprintf('%d (%1.1f%%)',sum(strcmp(lat_fmri,'L')),sum(strcmp(lat_fmri,'L'))/sum(fmri_non_control)*100)};
right_str = {'Right: N (%)',...
    sprintf('%d (%1.1f%%)',sum(right==1 & is_hup),sum(right==1 & is_hup)/sum(~isnan(right(is_hup)))*100),...
    sprintf('%d (%1.1f%%)',sum(right==1 & is_musc),sum(right==1 & is_musc)/sum(~isnan(right(is_musc)))*100),...
    sprintf('%d (%1.1f%%)',sum(strcmp(lat_fmri,'R')),sum(strcmp(lat_fmri,'R'))/sum(fmri_non_control)*100)};
bilat_str = {'Bilateral: N (%)',...
    sprintf('%d (%1.1f%%)',sum(bilateral==1 & is_hup),sum(bilateral==1 & is_hup)/sum((is_hup))*100),...
    sprintf('%d (%1.1f%%)',sum(bilateral==1 & is_musc),sum(bilateral==1 & is_musc)/sum((is_musc))*100),...
    sprintf('%d (%1.1f%%)',sum(strcmp(lat_fmri,'B')),sum(strcmp(lat_fmri,'B'))/sum(fmri_non_control)*100)};
surg_str = {'Surgery performed','','',''};
resection_str = {'Resection: N (%)',...
    sprintf('%d (%1.1f%%)',sum(resection==1 & is_hup),sum(resection==1 & is_hup)/sum(is_hup)*100),...
    sprintf('%d (%1.1f%%)',sum(resection==1 & is_musc),sum(resection==1 & is_musc)/sum(is_musc)*100),...
    sprintf('%d (%1.1f%%)',nresection_fmri,nresection_fmri/sum(fmri_non_control)*100)};
ablation_str = {'Ablation: N (%)',...
    sprintf('%d (%1.1f%%)',sum(ablation==1 & is_hup),sum(ablation==1 & is_hup)/sum(is_hup)*100),...
    sprintf('%d (%1.1f%%)',sum(ablation==1 & is_musc),sum(ablation==1 & is_musc)/sum(is_musc)*100),...
    sprintf('%d (%1.1f%%)',nablation_fmri,nablation_fmri/sum(fmri_non_control)*100)};
device_str = {'Device: N (%)',...
    sprintf('%d (%1.1f%%)',sum(device==1 & is_hup),sum(device==1 & is_hup)/sum(is_hup)*100),...
    sprintf('%d (%1.1f%%)',sum(device==1 & is_musc),sum(device==1 & is_musc)/sum(is_musc)*100),...
    sprintf('%d (%1.1f%%)',ndevice_fmri,ndevice_fmri/sum(fmri_non_control)*100)};
engel_str = {'Engel 1 year outcome',...
    sprintf('N = %d with outcomes',sum(~isnan(engel_yr1(is_hup&(resection | ablation))))),...
    sprintf('N = %d with outcomes',sum(~isnan(engel_yr1(is_musc&(resection | ablation))))),...
    sprintf('N = %d with outcomes',sum(~isnan(engel_yr1_fmri)))};
engel_one_str = {'Median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(engel_yr1(is_hup&(resection | ablation))),min(engel_yr1(is_hup&(resection | ablation))),max(engel_yr1(is_hup&(resection | ablation)))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(engel_yr1(is_musc)),min(engel_yr1(is_musc)),max(engel_yr1(is_musc))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(engel_yr1_fmri),min(engel_yr1_fmri),max(engel_yr1_fmri))};
engel_two_str = {'Year 2: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(engel_yr2(is_hup)),min(engel_yr2(is_hup)),max(engel_yr2(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(engel_yr2(is_musc)),min(engel_yr2(is_musc)),max(engel_yr2(is_musc))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(engel_yr2_fmri),min(engel_yr2_fmri),max(engel_yr2_fmri))};
ilae_str = {'ILAE 1 year outcome',...
    sprintf('N = %d with outcomes',sum(~isnan(ilae_yr1(is_hup&(resection | ablation))))),...
    sprintf('N = %d with outcomes',sum(~isnan(ilae_yr1(is_musc&(resection | ablation))))),...
    sprintf('N = %d with outcomes',sum(~isnan(ilae_yr1_fmri)))};
ilae_one_str = {'Median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(ilae_yr1(is_hup)),min(ilae_yr1(is_hup)),max(ilae_yr1(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(ilae_yr1(is_musc)),min(ilae_yr1(is_musc)),max(ilae_yr1(is_musc))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(ilae_yr1_fmri),min(ilae_yr1_fmri),max(ilae_yr1_fmri))};
ilae_two_str = {'Year 2: median (range)',...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(ilae_yr2(is_hup)),min(ilae_yr2(is_hup)),max(ilae_yr2(is_hup))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(ilae_yr2(is_musc)),min(ilae_yr2(is_musc)),max(ilae_yr2(is_musc))),...
    sprintf('%1.1f (%1.1f-%1.1f)',...
    nanmedian(ilae_yr2_fmri),min(ilae_yr2_fmri),max(ilae_yr2_fmri))};

all = [total_str;...
    female_str;...
    age_onset_str;...
    age_implant_str;...
    n_discordant_str;...
    n_elecs_str;...
   % loc_str;...
   % tle_str;...
   % etle_str;...
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
    %engel_two_str;...
    ilae_str;...
    ilae_one_str];
    %ilae_two_str];

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

%% Get spike detector performance
sT = readtable('Manual validation.xlsx','Sheet','EDF pipeline');
snames = sT.name;
bi_good = sT.x_Correct_outOf50__bi_;
perc_good = bi_good/50;

all_pts_sp = ismember(snames,T.names);
all_pts_hup_sp = ismember(snames,T.names) & contains(snames,'HUP');
all_pts_musc_sp = ismember(snames,T.names) & contains(snames,'MP');


%% Summary results
%fprintf(fid,'<p><br><u><i>Clinical information and summary of intracranial recording</i></u></br>');
%{
fprintf(fid,['We examined %d consecutive patients at HUP for our univariate analyses, model '...
    'development, and internal validation. We excluded %d patients without symmetric bitemporal '...
    'electrode coverage, and an additional %d patients with '...
    'extratemporal seizure onset zones, yielding %d patients included for analysis.'],...
    n_hup_total,n_exclude_no_symm_hup,n_etle_hup,n_hup_remaining);

fprintf(fid,[' We examined %d consecutive patients at MUSC for external model validation. We excluded %d patients '...
    'without symmetric bitemporal electrode coverage, an additional %d patient with extra-temporal seizure onset, '...
    'and an additional %d patient for whom most of the EEG was disconnected, yielding %d MUSC patients who '...
    'met inclusion criteria (Table 1).'],...
    n_musc_total,n_exclude_no_symm_musc,n_etle_musc,n_exclude_no_sleep_conn_musc,n_musc_remaining);
%}
fprintf(fid,['<p>We examined %d patients at HUP for model development and internal validation, '...
    'and %d patients at MUSC for external model validation (Table 1).'],n_hup_remaining,n_musc_remaining);
%{
fprintf(fid,[' The majority of patients underwent bilateral electrode implantation '...
    ' either because one of the two primary pre-implantation hypotheses was bilateral, '...
    'or because of discordant lateralities between the two primary pre-implant hypotheses.']);
%}

fprintf(fid,[' We visually validated a random sample of 50 automated spike detections from each patient '...
    '(bipolar montage). The median (IQR) percentage of automatically-detected spikes '...
    'determined to be true spikes was '...
    '%1.1f%% (%1.1f%%-%1.1f%%) for HUP and %1.1f%% (%1.1f%%-%1.1f%%) for MUSC.</p>'],...
    median(perc_good(all_pts_hup_sp))*100,prctile(perc_good(all_pts_hup_sp),25)*100,...
    prctile(perc_good(all_pts_hup_sp),75)*100,...
    median(perc_good(all_pts_musc_sp))*100,prctile(perc_good(all_pts_musc_sp),25)*100,...
    prctile(perc_good(all_pts_musc_sp),75)*100);

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

function engel = get_engel(tengel)

engel = cell(length(tengel),1);
for i = 1:length(tengel)
    switch tengel(i)
                    
        case 1
            engel{i} = 'IA';
        case 5
            engel{i} = 'IB';
        case 6
            engel{i} = 'IC';
        case 7
            engel{i} = 'ID';
        case 2
            engel{i} = 'IIA';
        case 8
            engel{i} = 'IIB';
        case 9
            engel{i} = 'IIC';
        case 10
            engel{i} = 'IID';
        case 3
            engel{i} = 'IIIA';
        case 11
            engel{i} = 'IIIB';
        case 4
            engel{i} = 'IVA';
        case 12
            engel{i} = 'IVB';
        case 13
            engel{i} = 'IVC';
        otherwise
            engel{i} = '';
            
    end

end
end

function ilae = get_ilae(tilae)
ilae = cell(length(tilae),1);
for i = 1:length(ilae)
    switch tilae(i)
                    
        case 1
            ilae{i} = 'ILAE 1';
        case 2
            ilae{i} = 'ILAE 1a';
        case 3
            ilae{i} = 'ILAE 2';
        case 4
            ilae{i} = 'ILAE 3';
        case 5
            ilae{i} = 'ILAE 4';
        case 6
            ilae{i} = 'ILAE 5';
        case 7
            ilae{i} = 'ILAE 6';
        otherwise
            ilae{i} = '';
    end
end

end