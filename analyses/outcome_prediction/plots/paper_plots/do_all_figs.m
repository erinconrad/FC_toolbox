function do_all_figs

%% Parameters
freqs = {'Delta (0.5-4 Hz)','Theta (4-8 Hz)','Alpha (8-12 Hz)','Beta (12-30 Hz)','Gamma (30-80 Hz)'};
def_colors = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.9290, 0.6940, 0.1250];

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];

bct_folder= locations.bct;
out_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/paper_plots/'];
model_folder = [results_folder,'analysis/outcome/plots/'];

if ~exist(plot_folder,'dir'), mkdir(plot_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));


%% Do all analyses
% Load results from model
model = load([model_folder,'model_stuff.mat']);

% symmetric coverage tests for both atlases
aal_out = symmetric_coverage_tests('aal_bernabei');
brain_out = symmetric_coverage_tests('brainnetome');

% Spike-fc corr
corr_out = spike_fc_correlation;



%% Figure 2 - symmetric coverage test for brainnetome and Supplemental Fig 1 (same but AAL)
% NEED TO THINK ABOUT N
% Prep section in text
for ia = 1:2
    
    if ia == 1
        atlas_out = brain_out;
        fid = fopen([plot_folder,'results.html'],'a');
        fprintf(fid,'<b>Connectivity with symmetric coverage constraint</b>');
        fig_name = 'Fig2';
    else
        atlas_out = aal_out;
        fid = fopen([plot_folder,'supplemental_results.html'],'a');
        fprintf(fid,'<b>Connectivity with symmetric coverage constraint - AAL atlas</b>');
        fig_name = 'FigS1';
    end
    
    

    fprintf(fid,['To control for spatial sampling bias, we first compared connectivity '...
        'in the SOZ against that in the contralateral region while restricting analysis '...
        'to those regions with symmetric electrode coverage.']);

    % unpack brain_out - be careful not to confuse with aal results
    unpack_any_struct(atlas_out);
    figure
    set(gcf,'position',[1 1 1400 900])
    tiledlayout(2,6,'tilespacing','tight','padding','tight')

    % Show bilateral coverage
    nexttile([1 3])
    pretty_matrix(all_bilateral(~neither_lat,:),...
        {'SOZ\newlinehemisphere','non-SOZ\newlinehemisphere'},sum(left),[],1);
    colormap(gca,[0.5,0.5,0.5;1 1 1]);
    title('Regions with symmetric coverage')
    xlabel('Patient')
    ylabel('Region')

    fprintf(fid,[' Across patients, there were on average %1.1f (range %d-%d) regions '...
        ' (%1.1f when separately counting left and right) with bilateral electrode coverage. '...
        '%d of %d patients had any regions with bilateral electrode coverage (Fig. 1A).'],...
        mean(sum(all_bilateral,1)/2),min(sum(all_bilateral,1)/2), max(sum(all_bilateral,1)/2),...
        mean(sum(all_bilateral,1)),sum(sum(all_bilateral,1)>0),size(all_bilateral,2));

    % why is n different?

    % Show atlas
    nexttile([1 3])
    pretty_matrix(nanmean(soz_non_soz_ordered_atlas(~neither_lat,~neither_lat,:),3),...
        {'SOZ\newlinehemisphere','non-SOZ\newlinehemisphere'},sum(left),'r',0);
    caxis(gca,[-1 1])
    title('Average connectivity (symmetric coverage only)')
    xlabel('Region')

    fprintf(fid,[' Visually, in a patient-averaged network, network edges in the hemisphere opposite the SOZ were '...
        'stronger than those in the hemisphere of the SOZ (Fig. 1B, brighter colors indicate higher connectivity).']);

    % Show tests
    %need TO ADD N FOR EACH TEST
    nexttile([1 2])
    stats = paired_plot(soz_all,'Average connectivity',{'SOZ','contralateral region','contralateral region'});
    title('Connectivity to SOZ vs. contralateral region')

    fprintf(fid,[' The average connectivity to the SOZ (median %1.1e) was lower than that to '...
        ' the region contralateral to the SOZ (median %1.1e) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 1C).'],...
        stats.medians(1),stats.medians(2),stats.Tpos,get_p_html(stats.pval));

    nexttile([1 2])
    stats = paired_plot(hemi,'Intra-hemispheric connectivity',{'SOZ side','non-SOZ side','non-SOZ side'});
    title({'Hemispheric connectivity','ipsilateral vs. contralateral to SOZ'})

    fprintf(fid,[' The average intra-hemispheric connectivity in the side of the SOZ (median %1.1e) was also lower than that in '...
        ' the contralateral hemisphere (median %1.1e) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 1D).'...
        ' This suggests that the change in connectivity in epilepsy is broad and affects the entire hemisphere.'],...
        stats.medians(1),stats.medians(2),stats.Tpos,get_p_html(stats.pval));

    nexttile([1 2])
    stats = paired_plot(soz_intra,'Intrinsic connectivity',{'SOZ','contralateral region','contralateral region'});
    title({'Intrinsic connectivity','in SOZ vs contralateral region'})

    fprintf(fid,[' In the case of patients whose SOZ spanned multiple regions, the intrinisic connectivity within the SOZ '...
        '(between one SOZ region and other SOZ regions) (median %1.1e) was no different from the intrinsic connectivity '...
        ' within the same regions in the contralateral hemisphere (median %1.1e) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 1E).'...
        ' This suggests no focal within-SOZ change in connectivity within the limits of this analysis.'],...
        stats.medians(1),stats.medians(2),stats.Tpos,get_p_html(stats.pval));

    % Add annotations
    annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
    annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',30,'linestyle','none')
    annotation('textbox',[0 0.44 0.1 0.1],'String','C','fontsize',30,'linestyle','none')
    annotation('textbox',[0.32 0.44 0.1 0.1],'String','D','fontsize',30,'linestyle','none')
    annotation('textbox',[0.66 0.44 0.1 0.1],'String','E','fontsize',30,'linestyle','none')

    print(gcf,[plot_folder,fig_name],'-dpng')
    
    fclose(fid);
end

%% Supplemental Fig 2 - coherence (both atlases)
figure
set(gcf,'position',[1 1 1440 900])
tiledlayout(2,5,'tilespacing','tight','padding','tight')
fig_name = 'FigS2';
fid = fopen([plot_folder,'supplemental_results.html'],'a');

for ia = 1:2
    
    if ia == 1
        atlas_out = brain_out;
        fprintf(fid,'<b>Connectivity with symmetric coverage constraint</b>');
        
    else
        atlas_out = aal_out;  
        fprintf(fid,'<b>Connectivity with symmetric coverage constraint - AAL atlas</b>');
    end
    
    % unpack brain_out - be careful not to confuse with aal results
    unpack_any_struct(atlas_out);
    
    for f = 1:5
        
        nexttile
        stats = paired_plot(squeeze(soz_coh_all(:,f,:)),'Coherence',{'SOZ','contralateral region','contralateral region'});
        title(freqs{f})

        fprintf(fid,[' The average %s coherence to the SOZ (median %1.1e) was lower than that to '...
            ' the region contralateral to the SOZ (median %1.1e) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 1C).'],...
            freqs{f},stats.medians(1),stats.medians(2),stats.Tpos,get_p_html(stats.pval));
        
    end

end

% Add annotations
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.2 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0.4 0.91 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.6 0.91 0.1 0.1],'String','D','fontsize',25,'linestyle','none')
annotation('textbox',[0.8 0.91 0.1 0.1],'String','E','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.41 0.1 0.1],'String','F','fontsize',25,'linestyle','none')
annotation('textbox',[0.2 0.41 0.1 0.1],'String','G','fontsize',25,'linestyle','none')
annotation('textbox',[0.4 0.41 0.1 0.1],'String','H','fontsize',25,'linestyle','none')
annotation('textbox',[0.6 0.41 0.1 0.1],'String','I','fontsize',25,'linestyle','none')
annotation('textbox',[0.8 0.41 0.1 0.1],'String','J','fontsize',25,'linestyle','none')
print(gcf,[plot_folder,fig_name],'-dpng')
fclose(fid);

%% Figure 3- confusion matrixes


for ia = 1:2
    
    figure
    set(gcf,'position',[1 100 800 370])
    tiledlayout(1,2,'tilespacing','compact','padding','tight');
    
    if ia == 1
        atlas_out = brain_out;
        fig_name = 'Fig3';
        
    else
        atlas_out = aal_out;  
        fig_name = 'FigS3';
    end
    
    unpack_any_struct(atlas_out);

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
    annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
    annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
    print(gcf,[plot_folder,fig_name],'-dpng')
end

%% Figure 4 - correlation between spikes and FC
unpack_any_struct(corr_out);

figure
set(gcf,'position',[1 100 800 370])
tiledlayout(1,2,'tilespacing','compact','padding','tight');

for i = 1:2
    if i == 1
        thing = spike_fc_corr;
        ttext = 'Raw connectivity';
    else
        thing = spike_resid_corr;
        ttext = 'Connectivity (normalized for density)';
    end
    nexttile
    pp = plot(thing,'o','linewidth',2);
    hold on
    plot(xlim,[0 0],'k--','linewidth',2)
    ylim([-1 1])
    median_corr = nanmedian(thing);
    xlabel('Patient')
    xticklabels([])
    ylabel('Spike rate - connectivity correlation (\rho)')
    title(ttext)
    [~,p] = ttest(thing);
    legend(pp,sprintf('median \\rho = %1.2f, %s',median_corr,get_p_text(p)),'location','southeast','fontsize',15)
    set(gca,'fontsize',15)
end
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
print(gcf,[plot_folder,'Fig4'],'-dpng')

%% Figure 5 - null model construction

%% Figure 6 - model
% Add conceptual figure
model_info = model.all_out.model_info;
mnames = {sprintf('Chance: AUC 0.50'),...
    sprintf('Null model: AUC %1.2f',mean(model_info(2).all_auc)),...
    sprintf('Connectivity: AUC %1.2f',mean(model_info(3).all_auc)),...
    sprintf('Spikes: AUC %1.2f',mean(model_info(4).all_auc))};
mbars = {'-','--','-.'};
mleg = nan(3,1);
figure
plot([0 1],[0 1],'k:','linewidth',2)
hold on
count = 0;
for im = [2 3 4]
    count = count+1;
    if im == 2
        mname = 'Null model';
    elseif im == 2
        mname = 'Connectivity';
    elseif im == 3
        mname = 'Spikes';
    end
    mleg(count) = plot(model_info(im).x,model_info(im).ym,'linestyle',mbars{count},'linewidth',2,...
        'color',def_colors(count,:));
    
    
end
set(gca,'fontsize',15)
xlabel('False positive rate')
ylabel('True positive rate')
legend(mnames,'fontsize',15,'location','southeast')
title('SOZ prediction')

end