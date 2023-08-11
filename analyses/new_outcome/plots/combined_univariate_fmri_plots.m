function combined_univariate_fmri_plots

%% Parameters
% Univariate parameters
which_pts = 'hup';
rm_non_temporal = 1;
response = 'soz_lats';
just_sep_bilat = 1;

% fmri parameters
rm_controls = 1;

%% Get file locs
locations = fc_toolbox_locs;
data_folder = locations.el_data_folder;
fmri_folder = [data_folder,'fmri_data/'];
plot_folder = locations.el_plots_folder;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load the file containing intermediate data
inter_folder = data_folder;
mt_data = load([inter_folder,'mt_out_epilepsy_laterality.mat']);
mt_data = mt_data.out;

%% Initialize results file
fname = [plot_folder,'results.html'];
fid = fopen(fname,'a');
fprintf(fid,'<br><u><i>Several interictal EEG feature asymmetries are distinct across SOZ lateralities</i></u></br>');
fprintf(fid,['We compared interictal EEG feature AIs '...
    'between patients with left-sided, right-sided, and bilateral SOZs.']);

%% Now do lr_mt to get AI features
[T,features] =  lr_mt(mt_data,3); % just sleep
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
xlim([0 n_to_plot+1])
ylim([min(feature_eta2(I(1:n_to_plot)))-0.02 max(feature_eta2(I(1:n_to_plot)))+0.02])
xticklabels(cellfun(all_shorten,cellfun(@greek_letters_plots,allowed_features(I(1:n_to_plot)),'uniformoutput',false),...
    'UniformOutput',false));
axt = gca;
axt.XAxisLocation = 'top';
axt.XTickLabelRotation = 45;
%}
title({'Effect size (\eta^2_p) to distinguish left/right/bilateral SOZ'})
set(gca,'fontsize',15)

fprintf(fid,[' We ranked '...
    'features in descending order by effect size (&#951;<sup>2</sup>) at '...
    'separating the three SOZ lateralities (Fig. 2A). '...
    'The top-ranked AI '...
    'features involved spike rates and '...
    'relative entropy. ']);

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
ax1.XLim = [0 n_to_plot+1];
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
ax2.XLim = [0 n_to_plot+1];

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
annotation('textbox', [0.55 0.895 0.1 0.1],...
    'String', 'Effect sizes (Cohen''s {\it d}) to distinguish specific laterality', ...
    'EdgeColor', 'none', ...
    'fontweight','bold','fontsize',17,...
    'HorizontalAlignment', 'left')

fprintf(fid,['We next compared the set of features that best distinguished '...
    'left from bilateral SOZs versus right from bilateral SOZs (Fig. 2B). '...
    'Spikes and relative entropy features significantly distinguished left from bilateral SOZs.'...
    ' For distinguishing '...
    'right-sided SOZ, several other features performed best (though none were significant'...
    ' after correcting for the false discovery rate). This suggests that right-sided SOZs '...
    'are harder to distinguish than left-sided SOZs in our dataset. '...
    'Several interictal features were highly correlated (see Fig. S1 and Supplemental Results).</p>']);


%% fmri locs
%file_path = [locations.main_folder,'Alfredo_code/fmri_analysis_AL_3_28_23/'];
file_path = fmri_folder;
csv_path = [file_path,'out_csvs/'];

%% Load files
T = readtable([file_path,'df.csv']);
bT = readtable([file_path,'BNA_subregions.xlsx']);
mt_mask = niftiread([file_path,'mesial_temporal_roi.nii.gz']);
mni_brain = niftiread([file_path,'tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz']);

fprintf(fid,'<br><u><i>fMRI BOLD connectivity AI also distinguishes SOZ lateralities</i></u></br>');
fprintf(fid,['We next asked if differences in interictal connectivity between TLE lateralities '...
    'existed across modalities beyond just EEG. We tested whether fMRI connectivity '...
    'asymmetry similarly distinguished left, right, '...
    'and bilateral SOZs.']);

%% Remove controls
if rm_controls
    controls = strcmp(T.Final_Lat,'Control');
    T(controls,:) = [];
end

%% Remove ieeg?
ieeg = strcmp(T.IEEG,'IEEG');


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

annotation('textbox', [0.07 0.32 0.1 0.1],...
    'String', 'Regions included in fMRI connectivity calculations', ...
    'EdgeColor', 'none', ...
    'fontweight','bold','fontsize',17,...
    'HorizontalAlignment', 'left')

fprintf(fid,[' Fig. 2C shows the left temporal Brainnetome atlas regions included for fMRI connectivity '...
    'analysis.']);


%% Main plot
nexttile(t,7,[1 2])
if rm_controls
    [~,stats] = boxplot_with_points(AI,lat,1,{'left','right','bilateral'},[],'para');
else
    [~,stats] = boxplot_with_points(AI,lat,1,{'left','right','bilateral','Control'},[],'para');
end
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

fprintf(fid,[' There was a significant difference in fMRI connectivity AI between TLE lateralities '...
    '(ANOVA: F(%d,%d) = %1.1f, %s, &#951;<sup>2</sup> = %1.2f, Fig. 2D). Only '...
    'the difference between the left and bilateral SOZ group was significant after '...
    'correcting for multiple comparisons (%s). '],...
    stats.tbl{2,3},stats.tbl{3,3},stats.tbl{2,5},get_p_html(stats.p),stats.eta2,...
    get_p_html(stats.lbp));

fprintf(fid,['Similar to the result for interictal EEG data, this suggests that fMRI temporal lobe connectivity '...
    'AI distinguishes patients with left from bilateral SOZs, but cannot '...
    'distinguish between patients with right and bilateral SOZs. ']);

% Get stats if remove ieeg
[p_rm_ieeg,tbl_rm_ieeg,stats_rm_ieeg] = anova1(AI(~ieeg),lat(~ieeg),'off');
eta2_rm_ieeg = tbl_rm_ieeg{2,2}/(tbl_rm_ieeg{2,2}+tbl_rm_ieeg{3,2});

%{
fprintf(fid,['Of the %d patients and controls studied in the fMRI analysis, '...
    '%d (%1.1f%%) also underwent intracranial EEG recording. If we exclude these patients '...
    'from our analysis, the difference in fMRI connectivity AI across lateralities is no longer significant'...
    ', although the effect size is similar (ANOVA: F(%d,%d) = %1.1f, %s, &#951;<sup>2</sup> = %1.2f).</p>'],...
    length(ieeg),sum(ieeg==1),sum(ieeg==1)/length(ieeg)*100,...
    tbl_rm_ieeg{2,3},tbl_rm_ieeg{3,3},tbl_rm_ieeg{2,5},get_p_html(p_rm_ieeg),eta2_rm_ieeg);
%}


%% Add subtitles
annotation('textbox',[0 0.9 0.1 0.1],'String','A','LineStyle','none','fontsize',25)
annotation('textbox',[0.49 0.9 0.1 0.1],'String','B','LineStyle','none','fontsize',25)
annotation('textbox',[0 0.33 0.1 0.1],'String','C','LineStyle','none','fontsize',25)
annotation('textbox',[0.49 0.33 0.1 0.1],'String','D','LineStyle','none','fontsize',25)


print(gcf,[plot_folder,'Fig2'],'-dpng')

%% Prediction analysis
%{
% Do LOO CV
npts = length(AI);

% Prepare combined lats
comb_br = lat;
comb_br(contains(comb_br,'bilateral')|contains(comb_br,'right')) = {'br'};
pos_br = 'left';

comb_bl = lat;
comb_bl(contains(comb_bl,'bilateral')|contains(comb_bl,'left')) = {'bl'};
pos_bl = 'right';

% Initialize scores
all_scores = nan(2,npts);
alt_scores = nan(2,npts);

for ic = 1:2

    % Choose which lats to combine
    if ic == 1
        comb_lat = comb_br;
        pos = pos_br;
    elseif ic == 2
        comb_lat = comb_bl;
        pos = pos_bl;
    end

    comb_lat_bin = zeros(npts,1);
    comb_lat_bin(strcmp(comb_lat,pos)) = 1;

    fT = table(comb_lat,comb_lat_bin,AI,'VariableNames',{'lat','lat_bin','AI'});

    % Loop over patients for LOO
    for ip = 1:npts

        % divide training and testing data
        testing = ip;
        training = [1:ip-1,ip+1:npts];

        ft_test = fT(testing,:);
        ft_train = fT(training,:);

        % train the model on the training data
        mdl = fitglm(ft_train,'lat_bin ~ AI','Distribution','binomial'); % basic LR model
        coef = [mdl.Coefficients.Estimate(1);mdl.Coefficients.Estimate(2)];

        % test the model on the testing data
        score = glmval(coef,ft_test.AI,'logit');
        alt_scores(ic,ip) = 1/(1+exp(-coef(1)-coef(2)*ft_test.AI));
        all_scores(ic,ip) = score;
       
        
    end
end

% Get ROC curves
[XL,YL,~,AUCL] = perfcurve(comb_br,all_scores(1,:),'left');
[XR,YR,~,AUCR] = perfcurve(comb_bl,all_scores(2,:),'bl');

% Plot ROC curves
figure
ll = plot(XL,YL,'linewidth',2);
hold on
%lr = plot(XR,YR,'linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
%{
legend([ll,lr],{sprintf('Left vs right/bilateral: AUC = %1.2f',AUCL),...
    sprintf('Right vs left/bilateral: AUC = %1.2f',AUCR)},'fontsize',15,...
    'location','southeast')
%}
legend([ll],{sprintf('Left vs right/bilateral: AUC = %1.2f',AUCL),...
    },'fontsize',15,...
    'location','southeast')
title({'LOO CV'})
%}


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

