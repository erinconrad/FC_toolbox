function fmri_plot

%{
Need to match lateralities with my lateralities
%}

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end
file_path = '/Users/erinconrad/Desktop/research/FC_toolbox/Alfredo_code/fmri_analysis_AL_3_28_23/';
csv_path = [file_path,'out_csvs/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load files
T = readtable([file_path,'df.csv']);
bT = readtable([file_path,'BNA_subregions.xlsx']);
mt_mask = niftiread([file_path,'mesial_temporal_roi.nii.gz']);
mni_brain = niftiread([file_path,'tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz']);

%% Define temporal ROIs
temporal_hippo_amygdala_left = [108 110 112 114 116 118  74  78  86 212 214 216];
temporal_hippo_amygdala_right = [109 111 113 115 117 119  75  79  87 213 215 217];

% Add one to the indices of the regions (python to matlab)
temporal_hippo_amygdala_left = temporal_hippo_amygdala_left + 1;
temporal_hippo_amygdala_right = temporal_hippo_amygdala_right + 1;

% Show the names of the regions
lids = bT.LabelID_L;
rids = bT.LabelID_R;
regions = bT.Var6;
alt_regions = bT.LeftAndRightHemisphere;

left_regions = alt_regions(ismember(lids,temporal_hippo_amygdala_left));
right_regions = alt_regions(ismember(rids,temporal_hippo_amygdala_right));
assert(isequal(left_regions,right_regions))

%% Get fcon
npatients = size(T,1);
all_fcon = nan(npatients,246,246);

% Loop over the patients
for i = 1:npatients
    % Get the subject id
    subj_id = T.Subject{i};

    % Load the csv containing the patient-specific fcon
    fcon_T = readtable([csv_path,subj_id,'.csv']);
    fcon = table2array(fcon_T);
    all_fcon(i,:,:) = fcon;

end

%% Get the left strength and right strength
left_str = mean(sum(abs(all_fcon(:,temporal_hippo_amygdala_left,temporal_hippo_amygdala_left)),3),2);
right_str = mean(sum(abs(all_fcon(:,temporal_hippo_amygdala_right,temporal_hippo_amygdala_right)),3),2);

%% Define AI
AI = (left_str-right_str)./(left_str+right_str);

%% Get lateralities
lat = T.Final_Lat;

%% Check lats
% Load manual validation file
mT = readtable('Manual validation.xlsx','Sheet','Comparing lats');


alats = cell(size(mT,1),1);
% loop over rows of this
for r = 1:size(mT,1)
    
    % grab rid
    rid = mT.RIDs(r);

    % turn it into format of T
    if rid < 100
        rid_text = sprintf('sub-RID00%d',rid);
    else
        rid_text = sprintf('sub-RID0%d',rid);
    end

    % find the matching row of Alfredo's table
    arow = strcmp(T.Subject,rid_text);

    if sum(arow) ~= 1, continue; end

    alats{r} = T.Final_Lat{arow};

end
alats(strcmp(alats,'R')) = {'right'};
alats(strcmp(alats,'L')) = {'left'};
alats(strcmp(alats,'B')) = {'bilateral'};

nmT = mT;
nmT = addvars(nmT,alats);

% remove empty
nmT(cellfun(@isempty,alats),:) = [];
assert(isequal(nmT.my_lats,nmT.alats))


% Now go the other way, what are the hup id numbers associated with
% Alfredo's RIDs?
hup_names = cell(size(T,1),1);
hup_lats = cell(size(T,1),1);
for r = 1:size(T,1)
    rid_text = T.Subject{r};
    % get just the number
    rid = str2num(strrep(rid_text,'sub-RID0',''));

    % find the matching RID in my lookup table
    mr = mT.RIDs == rid;

    if sum(mr) ~=1, continue; end
    hup_names{r} = mT.name{mr};
    hup_lats{r} = mT.my_lats{mr};
end

hup_lats(strcmp(hup_lats,'right')) = {'R'};
hup_lats(strcmp(hup_lats,'left')) = {'L'};
hup_lats(strcmp(hup_lats,'bilateral')) = {'B'};

naT = T;
naT = addvars(naT,hup_names);
naT = addvars(naT,hup_lats);

% remove empty
naT(cellfun(@isempty,hup_names),:) = [];
assert(isequal(naT.hup_lats,naT.Final_Lat))




% remove controls
controls = strcmp(lat,'Control');
lat(controls) = [];
AI(controls) = [];


% change names
lat(strcmp(lat,'L')) = {'left'};
lat(strcmp(lat,'R')) = {'right'};
lat(strcmp(lat,'B')) = {'bilateral'};

%% Make figure
figure
set(gcf,'position',[10 10 530 800])
tiledlayout(3,2,'tilespacing','tight','padding','compact')

%% Brains
special_color = [0.9 0.1 0.1];
% make it a double
mni_brain = double(mni_brain);

% resize mni brain to match the mt mask
mni_brain = imresize3(mni_brain,size(mt_mask));

% set 0 to be the brightest
mni_brain(mni_brain <= 100) = max(mni_brain,[],'all');

% Set the mt_mask to a special value
mni_brain(mt_mask==1) = nan;

% Plots
t1 = nexttile;
colormap(t1,'gray')
turn_nans_gray(imrotate(mni_brain(:,:,60),90),special_color,t1)
axis off
%title('Regions included in connectivity calculations')
set(gca,'fontsize',20)

% Plot brain again
t2 = nexttile;
colormap(t2,'gray')
turn_nans_gray(imrotate(squeeze(mni_brain(:,100,:)),90),special_color,t2)
axis off


%% Individual matrix plot
fcon1 = squeeze(all_fcon(30,:,:));
thing = [fcon1(temporal_hippo_amygdala_left,temporal_hippo_amygdala_left),fcon1(temporal_hippo_amygdala_left,temporal_hippo_amygdala_right);...
    fcon1(temporal_hippo_amygdala_right,temporal_hippo_amygdala_left),fcon1(temporal_hippo_amygdala_right,temporal_hippo_amygdala_right)];
nexttile([1 2])
turn_nans_gray(thing)
hold on
plot([size(thing,1)/2 size(thing,1)/2],[-2 size(thing,1)],'k-','linewidth',4)
plot([-1 size(thing,1)+0.5],[size(thing,1)/2 size(thing,1)/2],'k-','linewidth',4)
xticks([size(thing,1)/4 3*size(thing,1)/4])
yticks([size(thing,1)/4 3*size(thing,1)/4])
xticklabels({'Left','Right'})
yticklabels({'Left','Right'})
xlabel('Temporal region')
ylabel('Temporal region')
title('Single patient connectivity matrix')
set(gca,'fontsize',20)



%% Main plot
nexttile([1 2])
[~,stats] = boxplot_with_points(AI,lat,1,{'left','right','bilateral'});
set(gca,'fontsize',20)
ylabel('Asymmetry index')
title('fMRI asymmetry by SOZ laterality')

annotation('textbox', [0.05 0.90 1 0.1],...
    'String', 'Regions included in connectivity calculations', ...
    'EdgeColor', 'none', ...
    'fontweight','bold','fontsize',20,...
    'HorizontalAlignment', 'center')

print(gcf,[plot_folder,'Fig4'],'-dpng')
end