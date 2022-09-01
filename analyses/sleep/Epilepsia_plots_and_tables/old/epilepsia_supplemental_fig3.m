function epilepsia_supplemental_fig3

%% Parameters
myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];

%% Locations
locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/epilepsia/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

%% Prep supplementary results
fid = fopen([out_folder,'supplemental_results.html'],'a');
fprintf(fid,'<p><b>Model performance by implant coverage</b><br>');


%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;
npts = length(out.circ_out.names);
elec_locs = out.circ_out.all_elec_locs;
elec_lats = out.circ_out.all_elec_lats;
summ = count_locs_and_lats(elec_locs,elec_lats);

%% Get bilateral/unilateral and number of regions implanted
unilat = summ.loc_lat_count(:,2)==1;
n_non_empty = sum(summ.empty_locs_lats(:,1)==0);
bilat = summ.loc_lat_count(:,2)==2;
num_regions = summ.loc_lat_count(:,1);


%% Get model performances
nmout = out.model_out;
pt_stats = nmout.pt_stats;
pt_specific = nmout.pt_specific;
aucs = pt_stats(:,5);
scores =  pt_specific(:,1);
soz = pt_specific(:,3);
threshold = pt_specific(:,2);
desired_threshold = [];
[npv,ppv,pred,totaln,mat,threshold,acc] = cellfun(@(x,y,z) individual_threshold_stats(x,y,z,desired_threshold),...
    scores,soz,threshold,'uniformoutput',false);
ppv = cell2mat(ppv);

%% Initialize figure
figure
tiledlayout(2,1,'tilespacing','tight','padding','tight')

%% Compare model peformance based on bilateral/unilateral
nexttile
plot(1+randn(sum(unilat),1)*0.05,aucs(unilat),'o','linewidth',2,'color',myColours(1,:))
hold on
plot(2+randn(sum(bilat),1)*0.05,aucs(bilat),'o','linewidth',2,'color',[0.9290, 0.6940, 0.1250])
xticks([1 2])
xticklabels({'Unilateral','Bilateral'});
ylabel('Classification AUC')
yl = ylim;
ybar = yl(1) + 1.1*(yl(2)-yl(1));
ytext = yl(1) + 1.2*(yl(2)-yl(1));
nylim = [yl(1) yl(1) + 1.3*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k','linewidth',2)
text(1.5,ytext,get_p_text(ranksum(aucs(unilat),aucs(bilat))),...
    'fontsize',15,'horizontalalignment','center')
ylim(nylim)
xlim([0.5 2.5])
title('Classification performance by implant coverage')
set(gca,'fontsize',15)
[p,~,stats] = ranksum(aucs(unilat),aucs(bilat));
% text
W = stats.ranksum;
nt = sum(~(isnan(aucs(unilat))));
ne = sum(~(isnan(aucs(bilat))));
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);

fprintf(fid,['We tested whether the pre-implant hypothesis as measured by implant '...
    'coverage predicted model performance. We included patients with available '...
    'anatomical electrode localizations (N = %d). Most of the patients in the study '...
    'had bilateral electrode implantations with several regions sampled (Supplemental Table 1). '...
    'Model AUCs were similar between patients '...
    'with unilateral (median %1.2f) and bilateral (median %1.2f) implants '...
    '(Mann-Whitney test: <i>U</i>(<i>N<sub>unilateral</sub></i> = %d, <i>N<sub>bilateral</sub></i> = %d) ='...
    ' %1.1f, %s).'],n_non_empty,nanmedian(aucs(unilat)),nanmedian(aucs(bilat)),...
    nt,ne,U,get_p_html(p));



nexttile
plot(num_regions+randn(length(num_regions),1)*0.05,aucs,'o','linewidth',2)
hold on
[r,p] = corr(num_regions,aucs,'rows','pairwise','type','spearman');
yl = ylim;
ybar = yl(1) + 1.1*(yl(2)-yl(1));
ylabel('Classification AUC')
xlabel('Number of regions implanted')
ytext = yl(1) + 1.2*(yl(2)-yl(1));
nylim = [yl(1) yl(1) + 1.3*(yl(2)-yl(1))];
plot([1 6],[ybar ybar],'k','linewidth',2)
text(3.5,ytext,sprintf('\\rho = %1.2f, %s',r,get_p_text(p)),...
    'fontsize',15,'horizontalalignment','center')
ylim(nylim)
xlim([0.5 6.5])
set(gca,'fontsize',15)

fprintf(fid,[' Model AUCs were also not correlated with the number of regions implanted '...
    '(Spearman &rho; = %1.2f, %s) (Supplementary Fig. 3).'],r,get_p_html(p));

fclose(fid);
print(gcf,[out_folder,'FigS3'],'-depsc')
close all


end