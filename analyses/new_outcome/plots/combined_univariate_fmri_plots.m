function combined_univariate_fmri_plots

%% Parameters
% Univariate parameters
which_pts = 'hup';
rm_non_temporal = 1;
response = 'soz_lats';
just_sep_bilat = 1;
nboot = 1e1;

% fmri parameters
rm_ieeg = 0;
rm_controls = 0;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
subplot_path = [plot_folder,'ai_subplots/'];
if ~exist(subplot_path,'dir')
    mkdir(subplot_path)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Initialize results file
fname = [plot_folder,'results.html'];
fid = fopen(fname,'a');
fprintf(fid,'<p><br><b>Several interictal EEG feature asymmetries are distinct across SOZ lateralities</b></br>');
fprintf(fid,['We first examined, using univariate comparisons, how interictal EEG feature AI differs '...
    'between patients with left-sided SOZs, right-sided SOZs, and bilateral SOZs.']);

%% Now do lr_mt to get AI features
[T,features] =  lr_mt(3); % just sleep
allowed_features = features;

%% Restrict to desired hospital
switch which_pts
    case 'all'
    case 'hup'
        hup = contains(T.names,'HUP');
        T(~hup,:) = [];
    case 'musc'
        musc = contains(T.names,'MP');
        T(~musc,:) = [];
end

%% Remove non temporal patients
if rm_non_temporal
    temporal = strcmp(T.soz_locs,'temporal');
    T(~temporal,:) = [];
end

%% Initialize figure
figure
%set(gcf,'position',[25 235 1423 479])
%t = tiledlayout(2,4,"TileSpacing",'tight','padding','tight');

set(gcf,'position',[25 235 1050 900])
t = tiledlayout(2,4,"TileSpacing",'tight','padding','tight');


%% Univariate analyses of features
rm_sleep_text = @(x) strrep(x,' sleep','');
shorten_bi_text = @(x) strrep(x,'bipolar','bi');
shorten_machine_text = @(x) strrep(x,'machine','mac');
all_shorten = @(x) rm_sleep_text(shorten_machine_text(shorten_bi_text(x)));

nfeatures = length(allowed_features);
feature_p_val = nan(nfeatures,1);
feature_eta2= nan(nfeatures,1);

% Loop over features
for i = 1:nfeatures

    % get p value and eta 2
    [feature_p_val(i),tbl] = anova1(T.(allowed_features{i}),T.(response),'off');
    feature_eta2(i) = tbl{2,2}/(tbl{2,2}+tbl{3,2});
    
end
feature_p_val(isnan(feature_p_val)) = 1;

% False discovery rate
qvalues = mafdr(feature_p_val,'bhfdr',true); % Use BH approach, more conservative than Storey
n_to_plot = 15;
assert(sum(isnan(feature_eta2))==0)
[~,I] = sort(feature_eta2,'descend'); % sort by eta2
sorted_qvalues = qvalues(I(1:n_to_plot));

nexttile(t,1,[1 2])
plot(feature_eta2(I(1:n_to_plot)),'ko','markersize',15,'linewidth',2) % plot the top eta 2 scores
hold on
for i = 1:n_to_plot
    if sorted_qvalues(i) < 0.05
        plot(i,feature_eta2(I(i)),'k*','markersize',15,'linewidth',2)
    end
end
%assert(all(sorted_qvalues<0.05)) % ensure q is actually < 0.05

xticks(1:n_to_plot)
ylabel('\eta^2_p')

xticklabels(cellfun(all_shorten,cellfun(@greek_letters_plots,allowed_features(I(1:n_to_plot)),'uniformoutput',false),...
    'UniformOutput',false));
axt = gca;
axt.XAxisLocation = 'top';
axt.XTickLabelRotation = 45;
%}
title({'Effect size (\eta^2_p) to distinguish left/right/bilateral SOZ'})
set(gca,'fontsize',15)

fprintf(fid,['We examined the ability of interictal EEG features to '...
    'distinguish SOZ laterality. For each interictal EEG feature, we calculated the effect '...
    'size (&#951;<sup>2</sup>) at separating the three SOZ lateralities. We ranked '...
    'features in descending order by effect size. Fig. 3A shows the effect size of the top %d ranked features. '...
    'Each of these features had a significant effect at separating the three SOZ lateralities (ANOVA with '...
    'Benjamini-Hochberg false discovery rate correction). The top-ranked AI '...
    'features involve spike rates and '...
    'relative entropy.</p>'],n_to_plot);

%% Univariate analyses for different comparisons
feature_p_val = nan(nfeatures,2);
feature_eta2= nan(nfeatures,2);
for i = 1:nfeatures
    curr_feat = T.(allowed_features{i});

    % Just separate each from bilateral
    if just_sep_bilat
        [~,feature_p_val(i,1)] = ttest2(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'bilateral')));
        [~,feature_p_val(i,2)] = ttest2(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'bilateral')));
    
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,1) = dT.Effect;
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,2) = dT.Effect;
    else
        [~,feature_p_val(i,1)] = ttest2(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')));
        [~,feature_p_val(i,2)] = ttest2(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'left')|strcmp(T.(response),'bilateral')));
    
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,1) = dT.Effect;
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'left')|strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,2) = dT.Effect;
    end
end
%feature_p_val(isnan(feature_p_val)) = 1;
assert(sum(isnan(feature_eta2),'all')==0)
% False discovery rate
qvalues1 = mafdr(feature_p_val(:,1),'bhfdr',true); qvalues2 = mafdr(feature_p_val(:,2),'bhfdr',true);
n_to_plot = 15;
assert(sum(isnan(feature_eta2),'all')==0)
[~,I1] = sort(abs(feature_eta2(:,1)),'descend'); [~,I2] = sort(abs(feature_eta2(:,2)),'descend'); 
sorted_q1 = qvalues1(I1(1:n_to_plot)); sorted_q2 = qvalues2(I2(1:n_to_plot));
tt = tiledlayout(t,1,1,'tilespacing','none','padding','none');
tt.Layout.Tile = 3;
tt.Layout.TileSpan = [1 2];
ax1 = axes(tt);

pl = plot(ax1,1:n_to_plot,abs(feature_eta2(I1(1:n_to_plot),1)),'o','markersize',15,'color',[0 0.4470 0.7410],'linewidth',2);
hold on
for i =1:n_to_plot
    if sorted_q1(i) < 0.05
        plot(ax1,i,abs(feature_eta2(I1(i),1)),'*','markersize',15,'color',[0 0.4470 0.7410],'linewidth',2);
        
    end
end

ax1.XAxisLocation = 'top';
ax1.YAxisLocation = 'left';
ax1.XColor = [0 0.4470 0.7410];
ax1.YColor = [0 0.4470 0.7410];
ax1.Box = 'off';
ax1.YLim = [min([min(abs(feature_eta2(I1(1:n_to_plot),1))),min(abs(feature_eta2(I2(1:n_to_plot),2)))])-0.2,...
    max([max(abs(feature_eta2(I1(1:n_to_plot),1))),max(abs(feature_eta2(I2(1:n_to_plot),2)))])+0.2];
ax1.XTick = 1:n_to_plot;
%{
ax1.XTickLabel = cellfun(make_text_short,cellfun(@greek_letters_plots,allowed_features(I1(1:n_to_plot)),'uniformoutput',false),...
    'uniformoutput',false);
%}
%ax1.XTickLabel = cellfun(@make_text_short,allowed_features(I1(1:n_to_plot)),'uniformoutput',false);


ax2 = axes(tt);

pr = plot(ax2,abs(feature_eta2(I2(1:n_to_plot),2)),'o','markersize',15,'color',[0.8500 0.3250 0.0980],'linewidth',2);
hold on
for i =1:n_to_plot
    if sorted_q2(i) < 0.05
        plot(ax2,i,abs(feature_eta2(I2(i),2)),'*','markersize',15,'color',[0.8500 0.3250 0.0980],'linewidth',2);
    end
end

ax2.XAxisLocation = 'bottom';
ax2.YAxisLocation = 'right';
ax2.XColor = [0.8500 0.3250 0.0980];
ax2.YColor = [0.8500 0.3250 0.0980];

ax2.Color = 'none';
ax2.Box = 'off';

ax2.YLim = [min([min(abs(feature_eta2(I1(1:n_to_plot),1))),min(abs(feature_eta2(I2(1:n_to_plot),2)))])-0.2,...
    max([max(abs(feature_eta2(I1(1:n_to_plot),1))),max(abs(feature_eta2(I2(1:n_to_plot),2)))])+0.2];
if 0
    table(feature_eta2(I2,2),qvalues2(I2),feature_p_val(I2,2))
end

ax2.XTick = 1:n_to_plot;
%{
ax2.XTickLabel = cellfun(make_text_short,cellfun(@greek_letters_plots,allowed_features(I2(1:n_to_plot)),'uniformoutput',false),...
    'uniformoutput',false);
%}
%ax2.XTickLabel = cellfun(@make_text_short,allowed_features(I2(1:n_to_plot)),'uniformoutput',false);
ax2.XTickLabel = cellfun(all_shorten, ...
    cellfun(@greek_letters_plots,allowed_features(I2(1:n_to_plot)),'uniformoutput',false),...
    'uniformoutput',false);
ax1.XTickLabel = cellfun(all_shorten, ...
    cellfun(@greek_letters_plots,allowed_features(I1(1:n_to_plot)),'uniformoutput',false),...
    'uniformoutput',false);
ax1.XTickLabelRotation =45;
ax2.XTickLabelRotation =45;

%xticklabels([])
if just_sep_bilat
    legend([pl pr],{'Left vs bilateral','Right vs bilateral'},'location','northeast','fontsize',15)
else
    legend([pl pr],{'Left vs right/bilateral','Right vs left/bilateral'},'location','northeast','fontsize',15)
end

set(ax1,'fontsize',15); set(ax2,'fontsize',15)
ylabel(ax1,'|Cohen''s {\it d}|','color','k','fontsize',15)
%{
title(tt,{'Effect sizes (Cohen''s {\it d}) to distinguish specific laterality'},...
    'fontsize',20,'fontweight','bold')
%}
annotation('textbox', [0.52 0.9 0.1 0.1],...
    'String', 'Effect sizes (Cohen''s {\it d}) to distinguish specific laterality', ...
    'EdgeColor', 'none', ...
    'fontweight','bold','fontsize',17,...
    'HorizontalAlignment', 'left')

fprintf(fid,['<p>We compared the set of interictal features that best distinguished '...
    'left from bilateral SOZs versus right from bilateral SOZs. '...
    'We separately calculated the absolute value of the effect size (Cohen''s <i>d</i>) '...
    'at distinguishing left-sided SOZs from bilateral SOZs and that for distinguishing '...
    'right-sided SOZs from bilateral SOZs. For each of the two classification questions, '...
    'we ranked features in descending order by their absolute Cohen''s <i>d</i>. Fig. 3B '...
    'shows the Cohen''s <i>d</i> values for the top %d ranked features. '...
    'For distinguishing '...
    'left-sided SOZ, spikes and relative entropy features performed best '...
    '. For distinguishing '...
    'right-sided SOZ, several other features performed best (though none were significant'...
    ' correcting for the false discovery rate). This suggests that right-sided SOZs '...
    'are harder to distinguish than left-sided SOZs in our dataset.</p>'],n_to_plot);


%% fmri locs
file_path = '/Users/erinconrad/Desktop/research/FC_toolbox/Alfredo_code/fmri_analysis_AL_3_28_23/';
csv_path = [file_path,'out_csvs/'];

%% Load files
T = readtable([file_path,'df.csv']);
bT = readtable([file_path,'BNA_subregions.xlsx']);
mt_mask = niftiread([file_path,'mesial_temporal_roi.nii.gz']);
mni_brain = niftiread([file_path,'tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz']);

fprintf(fid,'<p><br><b>fMRI connectivity AI also distinguishes SOZ lateralities</b></br>');
fprintf(fid,['We asked whether the AI of fMRI connectivity similarly distinguished left, right, '...
    'and bilateral SOZs.']);

%% Remove controls
if rm_controls
    controls = strcmp(T.Final_Lat,'Control');
    T(controls,:) = [];
end

%% Remove ieeg?
ieeg = strcmp(T.IEEG,'IEEG');
if rm_ieeg
    
    T(ieeg,:) = [];
end


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


%% change names
lat(strcmp(lat,'L')) = {'left'};
lat(strcmp(lat,'R')) = {'right'};
lat(strcmp(lat,'B')) = {'bilateral'};

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
t1 = nexttile(t,5);
colormap(t1,'gray')
turn_nans_gray(imrotate(mni_brain(:,:,60),90),special_color,t1)
%axis equal
axis off
%title({'Regions included in fMRI','connectivity calculations'})
set(gca,'fontsize',15)

% Plot brain again
t2 = nexttile(t,6);
colormap(t2,'gray')
turn_nans_gray(imrotate(squeeze(mni_brain(:,100,:)),90),special_color,t2)
%axis equal
axis off

annotation('textbox', [0.07 0.30 0.1 0.1],...
    'String', 'Regions included in fMRI connectivity calculations', ...
    'EdgeColor', 'none', ...
    'fontweight','bold','fontsize',17,...
    'HorizontalAlignment', 'left')

fprintf(fid,[' We parcellated brain regions according to the DKT atlas, and we identified'...
    ' regions in the temporal lobe gray matter (we excluded white matter regions as '...
    'the DKT atlas does not specify the anatomical location of different white matter regions. '...
    'We measured the fMRI BOLD correlation between all ipsilateral temporal regions, separately for '...
    'the left and the right. The goal was to approximate our measure of intra-hemispheric '...
    'temporal lobe correlations from the interictal EEG data as closely as possible '...
    'using fMRI data. Fig. 3C shows the DKT atlas regions included for fMRI connectivity '...
    'analysis (only the left-sided regions are highlighted).']);


%% Main plot
nexttile(t,7,[1 2])
[~,stats] = boxplot_with_points(AI,lat,1,{'left','right','bilateral','Control'},[],'para');
set(gca,'fontsize',15)
ylabel('Asymmetry index')
title({'fMRI asymmetry by SOZ laterality'})

%{
annotation('textbox', [0.05 0.90 1 0.1],...
    'String', 'Regions included in connectivity calculations', ...
    'EdgeColor', 'none', ...
    'fontweight','bold','fontsize',20,...
    'HorizontalAlignment', 'center')
%}

fprintf(fid,[' We defined the fMRI connectivity AI using the same method as for '...
    'interictal EEG data: <i>AI</i> = (<i>Connectivity</i><sub>Left</sub> - '...
    '<i>Connectivity</i><sub>Right</sub>)/(<i>Connectivity</i><sub>Left</sub> + '...
    '<i>Connectivity</i><sub>Right</sub>). We compared the AI between '...
    'patients with left-sided SOZs, right-sided SOZs, bilateral SOZs, and control '...
    'subjects without epilepsy (Fig. 3D). There was a significant difference in AI between groups '...
    '(ANOVA: F(%d,%d) = %1.1f, %s, &#951;<sup>2</sup> = %1.2f). In post-hoc t-tests, only '...
    'the difference between the left and bilateral SOZ group was significant after '...
    'correcting for multiple comparisons using Bonferroni''s method (%s). This result suggests that, similar to the result for interictal EEG data, fMRI temporal lobe connectivity '...
    'AI can distinguish patients with left from bilateral SOZs, but cannot clearly '...
    'distinguish between patients with right and bilateral SOZs.</p>'],...
    stats.tbl{2,3},stats.tbl{3,3},stats.tbl{2,5},get_p_html(stats.p),stats.eta2,...
    get_p_html(stats.lbp));




%% Add subtitles
annotation('textbox',[0 0.9 0.1 0.1],'String','A','LineStyle','none','fontsize',20)
annotation('textbox',[0.5 0.9 0.1 0.1],'String','B','LineStyle','none','fontsize',20)
annotation('textbox',[0 0.31 0.1 0.1],'String','C','LineStyle','none','fontsize',20)
annotation('textbox',[0.5 0.31 0.1 0.1],'String','D','LineStyle','none','fontsize',20)


print(gcf,[plot_folder,'Fig3'],'-dpng')


end

function xout = make_text_short(x)

if contains(x,'spikes')
    xout = 'S';
elseif contains(x,'re')
    xout = 'R';
else
    xout = 'O';
end

end

