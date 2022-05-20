function do_all_figs(doing_from_github)

close all

%% Parameters
freqs = {'Delta (0.5-4 Hz)','Theta (4-8 Hz)','Alpha (8-12 Hz)','Beta (12-30 Hz)','Gamma (30-80 Hz)'};
def_colors = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.9290, 0.6940, 0.1250;...
    0.4940, 0.1840, 0.5560];

%% Get file locs
locations = fc_toolbox_locs;
plot_folder = locations.paper_plot_folder;

if ~exist(plot_folder,'dir'), mkdir(plot_folder); end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
model_folder = locations.paper_plot_folder;

%% Delete the current results file
if exist([plot_folder,'results.html'],'file')~=0
    delete([plot_folder,'results.html']);
end

if exist([plot_folder,'supplemental_results.html'],'file')~=0
    delete([plot_folder,'supplemental_results.html']);
end

%% Load intermediate datasets
% Load results from model
brain_model = load([model_folder,'model_stuff_brainnetome.mat']);
aal_model = load([model_folder,'model_stuff_aal_bernabei.mat']);

% symmetric coverage tests for both atlases
aal_out = load([model_folder,'symm_cov_aal_bernabei.mat']);
aal_out = aal_out.nout;
brain_out = load([model_folder,'symm_cov_brainnetome.mat']);
brain_out = brain_out.nout;

% Spike-fc corr
corr_out = load([model_folder,'spike_analysis.mat']);
corr_out = corr_out.nout;

% density out
dens_out = load([model_folder,'dens_model.mat']);
dens_out = dens_out.out;

%% Prep some general results
fid = fopen([plot_folder,'results.html'],'a');
fprintf(fid,['<p>We included all patients who had available electrode localizations (%d patients), although the number '...
    'of patients analyzed varied by analysis, as noted in the results of individual '...
    'analyses. Patients were heterogeneous by age, sex, seizure localization ',...
    'and lateralization, and implant strategy (Table 1).</p>'],sum(corr_out.pts_with_any_locs));
fclose(fid);

%% Fig 1 - conceptual fig
if doing_from_github == 0
    main_conceptual_figure
end

%% Supplemental fig 1 - density model construction
all_fc = dens_out.all_fc;
all_locs = dens_out.all_locs;
vec_dens = dens_out.vec_dens;
vec_conn = dens_out.vec_conn;
sr = dens_out.sr;
all_r2 = dens_out.all_r2;
poss_sr = dens_out.poss_sr;
g = dens_out.g;

    
figure
set(gcf,'position',[10 10 1000 700])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

nexttile
D = make_interdist_matrix(all_locs);
dens = interdistance_to_density_matrix(D,sr);
turn_nans_gray(dens)
c = colorbar;
ylabel(c,'Density (mm)')
set(gca,'fontsize',15)
xticklabels([])
yticklabels([])
xlabel('Electrode')
ylabel('Electrode')
title({'Inter-electrode density','(single patient)'})


% Functional connectivity
nexttile
conn = all_fc;
turn_nans_gray(conn)
c = colorbar;
caxis([-1 1])
ylabel(c,'Pearson correlation')
set(gca,'fontsize',15)
xticklabels([])
yticklabels([])
xlabel('Electrode')
ylabel('Electrode')
title({'Functional connectivity','(single patient)'})

nexttile
%[xbins,ybins] = bin_data_get_means(vec_dens,vec_conn,1e3);
plot(vec_dens,vec_conn,'o')
hold on
temp_x = [min(vec_dens):0.1*(max(vec_dens)-min(vec_dens)):max(vec_dens)];
%temp_y = (f.p1 * temp_x + f.p2)./(temp_x + f.q1);
temp_y = (g.p1 + temp_x*g.p2);
rp = plot(temp_x,temp_y,'k','linewidth',3);
ylim([-1 1])
xlim([min(vec_dens) max(vec_dens)])
xlabel('Density (mm)')
ylabel('Functional connectivity (r)')
set(gca,'fontsize',15)
legend(rp,sprintf('r^2 = %1.2f',max(all_r2)),'fontsize',15,'location','northwest')
title({'Density-functional connectivity correlation','(all patients)'})

nexttile
plot(poss_sr,all_r2,'linewidth',2)
hold on
rp = plot(sr,all_r2(poss_sr == sr),'*','markersize',20,'linewidth',2);
ylim([0 0.35])
xlabel('Search radius')
ylabel('Model r^2')
title({'Density-functional connectivity model r^2','as a function of search radius'})
set(gca,'fontsize',15)
legend(rp,'Optimal search radius','fontsize',15)


annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.42 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.42 0.1 0.1],'String','D','fontsize',25,'linestyle','none')

print(gcf,[plot_folder,'Fig S1'],'-dpng')

%% Figure 2 - symmetric coverage test for brainnetome and Supplemental Fig 2 (same but AAL)
% NEED TO THINK ABOUT N
% Prep section in text
for ia = 1:2
    
    if ia == 1
        atlas_out = brain_out;
        fid = fopen([plot_folder,'results.html'],'a');
        fprintf(fid,'<p><br><b>Connectivity with symmetric coverage constraint</b></br>');
        fig_name = 'Fig 2';
        
        fprintf(fid,['To control for spatial sampling bias, we first compared connectivity '...
        'in the SOZ against that in the contralateral region while restricting analysis '...
        'to regions with symmetric electrode coverage.']);
        
    else
        atlas_out = aal_out;
        fid = fopen([plot_folder,'supplemental_results.html'],'a');
        fprintf(fid,'<p><br><b>Connectivity with symmetric coverage constraint - AAL atlas</b></br>');
        fig_name = 'Fig S2';
        
        fprintf(fid,['We repeated our symmetric coverage analysis '...
            'using the AAL atlas rather than the Brainnetome atlas to '...
            'parcellate brain regions.']);
    end
    
    

    

    % unpack brain_out - be careful not to confuse with aal results
    %unpack_any_struct(atlas_out);
    figure
    set(gcf,'position',[1 1 1400 900])
    tiledlayout(2,6,'tilespacing','tight','padding','tight')

    % Show bilateral coverage
    nexttile([1 3])
    pretty_matrix(atlas_out.all_bilateral(~atlas_out.neither_lat,:),...
        {'SOZ\newlinehemisphere','non-SOZ\newlinehemisphere'},sum(atlas_out.left),[],1);
    colormap(gca,[0.5,0.5,0.5;1 1 1]);
    title('Regions with symmetric coverage')
    xlabel('Patient')
    ylabel('Region')

    fprintf(fid,[' Across patients, there were on average %1.1f (range %d-%d) regions '...
        ' (%1.1f when separately counting left and right) with bilateral electrode coverage. '...
        '%d of %d patients had any regions with bilateral electrode coverage (%sA).'],...
        mean(atlas_out.n_coverage.nsymmetric),min(atlas_out.n_coverage.nsymmetric),max(atlas_out.n_coverage.nsymmetric),...
        mean(atlas_out.n_coverage.nsymmetric)*2,atlas_out.n_coverage.any_symmetric,atlas_out.total_n_for_symmetric,fig_name);

    % Show atlas
    nexttile([1 3])
    pretty_matrix(nanmean(atlas_out.soz_non_soz_ordered_atlas(~atlas_out.neither_lat,~atlas_out.neither_lat,:),3),...
        {'SOZ\newlinehemisphere','non-SOZ\newlinehemisphere'},sum(atlas_out.left),'r',0);
    caxis(gca,[-1 1])
    title('Average connectivity (symmetric coverage only)')
    xlabel('Region')

    fprintf(fid,[' Visually, in a network averaged across patients, network edges in the hemisphere opposite the SOZ were '...
        'stronger than those in the hemisphere of the SOZ (%sB, brighter colors indicate higher connectivity).'],fig_name);

    % Show tests
    %need TO ADD N FOR EACH TEST
    nexttile([1 2])
    stats = paired_plot(atlas_out.soz_all,'Average connectivity',{'to SOZ','to contralateral region'});
    title('Connectivity to SOZ vs. contralateral region')
    
    fprintf(fid,[' We compared the average connectivity to the SOZ with that to the '...
        'region contralateral to the SOZ. We included patients who had electrode coverage of both the SOZ '...
        'and the contralateral region in this analysis (%d patients).'],...
        atlas_out.n_analyses.n_soz_all);

    fprintf(fid,[' The average connectivity to the SOZ (median %s) was lower than that to '...
        ' the region contralateral to the SOZ (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sC).'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name);

    nexttile([1 2])
    stats = paired_plot(atlas_out.hemi,'Intra-hemispheric connectivity',{'in SOZ side','in non-SOZ side'});
    title({'Hemispheric connectivity','ipsilateral vs. contralateral to SOZ'})
    
    fprintf(fid,[' We next compared the intra-hemispheric connectivity between '...
        'the side of the SOZ and the contralateral hemisphere. For this analysis we included '...
        'patients with unilateral epilepsy and at least two symmetrically-implanted atlas regions in '...
        'each hemisphere, which allowed us to calculate intra-hemispheric connectivity (%d patients).'],...
        atlas_out.n_analyses.n_holo_hem);

    fprintf(fid,[' The average intra-hemispheric connectivity on the side of the SOZ (median %s) was also lower than that in '...
        ' the contralateral hemisphere (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sD).'...
        ' This suggests that the change in connectivity in epilepsy is broad and affects the entire hemisphere.'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name);

    nexttile([1 2])
    stats = paired_plot(atlas_out.soz_intra,'Intrinsic connectivity',{'in SOZ','in contralateral region'});
    title({'Intrinsic connectivity','in SOZ vs contralateral region'})
    
    fprintf(fid,[' We next compared the intrinsic connectivity within the SOZ '...
        'and that of the contralateral region. For this analysis we included '...
        'patients whose SOZ spanned at least two atlas regions, and who had electrode coverage '...
        'of these regions and of the contralateral regions (%d patients).'],...
        atlas_out.n_soz_intra);

    fprintf(fid,[' The intrinisic connectivity within the SOZ '...
        '(between one SOZ region and other SOZ regions) (median %s) was no different from the intrinsic connectivity '...
        ' within the same regions in the contralateral hemisphere (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sE).'...
        ' This suggests no focal within-SOZ change in connectivity within the limits of this analysis.'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name);

    % Add annotations
    annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
    annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',30,'linestyle','none')
    annotation('textbox',[0 0.44 0.1 0.1],'String','C','fontsize',30,'linestyle','none')
    annotation('textbox',[0.32 0.44 0.1 0.1],'String','D','fontsize',30,'linestyle','none')
    annotation('textbox',[0.66 0.44 0.1 0.1],'String','E','fontsize',30,'linestyle','none')

    print(gcf,[plot_folder,fig_name],'-dpng')
    
    if ia == 1
        fprintf(fid,[' Results were similar when using the AAL rather than the Brainnetome '...
            'atlas for parcellating brain regions (Supplemental Results; %s).'],fig_name);
    end
    
    fprintf(fid,[' Results when studying coherence rather than Pearson correlation networks '...
        'were more heterogeneous, seen for different frequency bands in the two different atlases '...
        '(Supplemental Results; %s).</p>'],fig_name);
    fclose(fid);
end

%% Supplemental Fig 3 - coherence (both atlases)

figure
set(gcf,'position',[1 1 1440 900])
tiledlayout(2,5,'tilespacing','tight','padding','tight')
fig_name = 'Fig S3';
fid = fopen([plot_folder,'supplemental_results.html'],'a');
fprintf(fid,'<br><b>Coherence-based connectivity with symmetric coverage constraint</b></br>');
for ia = 1:2
    
    if ia == 1
        atlas_out = brain_out;
        fprintf(fid,'<p><br><b>Brainnetome atlas</b></br>');
        
    else
        atlas_out = aal_out;  
        fprintf(fid,'<br><b>AAL atlas</b></br>');
    end
    
    % unpack brain_out - be careful not to confuse with aal results
    unpack_any_struct(atlas_out);
    
    for f = 1:5
        
        nexttile
        if ia == 1 && f == 1
            stats = paired_plot(squeeze(soz_coh_all(:,f,:)),'Coherence',...
                {'SOZ','in contralateral region','in contralateral region'},...
                0,[0.045,0.85,0.16,0.08]);

        else
            stats = paired_plot(squeeze(soz_coh_all(:,f,:)),...
                'Coherence',{'SOZ','in contralateral region','in contralateral region'},...
                1);
        end
        if ia == 1
            title(sprintf('Brainnetome: %s',freqs{f}))
        else
            title(sprintf('AAL: %s',freqs{f}))
        end
        if stats.pval < 0.05
            fprintf(fid,[' The average %s coherence to the SOZ (median %s) was lower than that to '...
                ' the region contralateral to the SOZ (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sC).'],...
                freqs{f},pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name);
        else
            fprintf(fid,[' There was no significant difference between the %s coherence to the SOZ (median %s) and that to '...
                ' the region contralateral to the SOZ (median %1.1e) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sC).'],...
                freqs{f},pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name);
        end
        
    end
    
    fprintf(fid,'</p>');

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

%% Figure 3 and supplemental fig 4- confusion matrixes

for ia = 1:2
    
    figure
    set(gcf,'position',[1 100 800 370])
    tiledlayout(1,2,'tilespacing','compact','padding','tight');
    
    if ia == 1
        atlas_out = brain_out;
        fig_name = 'Fig 3';
        fid = fopen([plot_folder,'results.html'],'a');
        fprintf(fid,'<p><br><b>Epilepsy lateralization</b></br>');
        fprintf(fid,['We asked how well intra-hemispheric functional connectivity '...
            'could lateralize epilepsy.']);
    else
        atlas_out = aal_out;  
        fig_name = 'Fig S4';
        fid = fopen([plot_folder,'supplemental_results.html'],'a');
        fprintf(fid,'<p><br><b>Epilepsy lateralization - AAL atlas</b></br>');
        fprintf(fid,['We repeated the epilepsy lateralization analysis using the AAL '...
            'atlas rather than the Brainnetome atlas.']);
    end
    
    unpack_any_struct(atlas_out);
    
    fprintf(fid,[' We restricted analysis to patients with unilateral epilepsy '...
        'and at least two symmetrically-implanted atlas regions in '...
        'each hemisphere, which allowed us to calculate intra-hemispheric connectivity (%d patients).'],...
        atlas_out.n_analyses.n_holo_hem);
    
    fprintf(fid,[' We predicted that the hemisphere with lower average intrinsic connectivity '...
        'was the side of the SOZ.']);

    for i = 1:2 % connectivity then spikes

        if i == 1
            conf_out = conf_out_fc;
            ttext = 'Connectivity';
            fprintf(fid,[' Our prediction was %1.1f%% accurate (PPV for identifying right-sided epilepsy %1.1f%%, NPV %1.1f%%) '...
                'at lateralizing the SOZ.'],conf_out.accuracy*100,...
                conf_out.ppv*100,conf_out.npv*100);
        else
            conf_out = conf_out_spikes;
            ttext = 'Spikes';
            if ia == 1
                fprintf(fid,[' For comparison, a prediction using spike rates '...
                    '(predicting the hemisphere of the SOZ to be the hemisphere with more spikes, restricting regions to '...
                    'those with symmetric coverage) was %1.1f%% accurate (PPV for identifying right-sided epilepsy %1.1f%%, NPV %1.1f%%) '...
                    'at lateralizing the SOZ (%s).'],conf_out.accuracy*100,...
                    conf_out.ppv*100,conf_out.npv*100,fig_name);
            else
                fprintf(fid,[' For comparison, a prediction using spike rates '...
                    '(predicting the hemisphere of the SOZ to be the hemisphere with more spikes, restricting regions to '...
                    'those with symmetric coverage) was %1.1f%% accurate (PPV for identifying right-sided epilepsy %1.1f%%, NPV %1.1f%%) '...
                    'at lateralizing the SOZ (%s).'],conf_out.accuracy*100,...
                    conf_out.ppv*100,conf_out.npv*100,fig_name);
            end
            
            if ia == 1
                fprintf(fid,[' Results were similar using the AAL atlas '...
                    '(Supplemental Results, Fig S4).']);
            end
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
    
    fprintf(fid,[' These results indicate that reduced intra-hemispheric connectivity '...
        'can lateralize epilepsy.</p>']);
    fclose(fid);
end

%% Figure 4 - correlation between spikes and FC
unpack_any_struct(corr_out);

figure
set(gcf,'position',[1 100 1300 370])
tiledlayout(1,3,'tilespacing','compact','padding','tight');
fid = fopen([plot_folder,'results.html'],'a');
fig_name = 'Fig 4';
fprintf(fid,'<p><br><b>Correlation between functional connectivity and spikes</b></br>');

fprintf(fid,['We next measured the correlation between functional connecitivity and spike rates. '...
    'For each patient, we took the Spearman rank correlation between the average spike rate '...
    'for each electrode and the average functional connectivity of each electrode.']);
for i = 1
    if i == 1
        thing = spike_fc_corr;
        ttext = 'Raw connectivity';
    else
        thing = spike_resid_corr;
        ttext = 'Connectivity (normalized for density)';
    end
    nexttile
    thing(isnan(thing)) = [];
    assert(length(thing) == sum(good_spikes))
    fprintf(fid,['% For this analysis, we included all patients with accurate '...
        'spike detections (%d patients).'],length(thing));
    pp = plot(thing,'o','linewidth',2);
    hold on
    plot(xlim,[0 0],'k--','linewidth',2)
    xlim([0 sum(good_spikes)])
    ylim([-1 1])
    mean_corr = nanmean(thing);
    xlabel('Patient')
    xticklabels([])
    ylabel('Spike rate - connectivity correlation (\rho)')
    title('Correlation between spikes and connectivity')
    [~,p,~,stats] = ttest(thing);
    legend(pp,sprintf('mean \\rho = %1.2f, %s',mean_corr,get_p_text(p)),'location','southeast','fontsize',15)
    set(gca,'fontsize',15)
    
    fprintf(fid,[' The mean correlation was &rho; = %1.2f, which was significantly less than zero '...
    '(<i>t</i>-test, <i>t</i>(%d) = %1.2f, %s (%sA).'],mean_corr,stats.df,stats.tstat,get_p_html(p),fig_name);
    fprintf(fid,[' This implies that electrodes with higher spike rates tend to be less connected '...
        'to other electrodes.']);
    
    nexttile
    stats = paired_plot(all_all_chs_corr,'Average connectivity',{'before spike','during spike'});
    set(gca,'fontsize',15)
    title({'Change in connectivity with single spikes',' (All electrodes)'})
    
   fprintf(fid,[' In an attempt to understand if the negative relationship between spikes and connectivity '...
        'is due to the spikes themselves reducing connectivity, we compared the functional '...
        'connectivity during individual spikes (2 s window centered around each spike) to that before the spike '...
        '(the prior 2 s window). We examined 100 random spikes for each patient, averaging the pre- and during-spike '...
        'connectivity across all spikes for each patient. The average connectivity across all electrodes '...
        'during the spike was higher (median %s) than that before the spike (median %s) '...
        '(Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sB). This implies that average global network connectivity '...
        'is increased during a spike.'],pretty_exp_html(stats.medians(2)),pretty_exp_html(stats.medians(1)),stats.Tpos,...
        get_p_html(stats.pval),fig_name);
    
    nexttile
    stats = paired_plot(all_single_chs_corr,'Average connectivity',{'before spike','during spike'});
    title({'Change in connectivity with single spikes',' (Spike electrode)'})
    set(gca,'fontsize',15)
    
    fprintf(fid,[' In order to examine the local effect of spikes, we also compared the average pre- and during-spike '...
        'for the single electrode the spike was detected on. The average connectivity in the spike electrode '...
        'during the spike was lower (median %s) than that before the spike (median %s) '...
        '(Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sB). '...
        'This implies that the average connectivity to the spike electrode is reduced during a spike. '...
        'The reduction in connectivity in the spike electrode relative to other electrodes '...
        'also suggests that the negative correlation between spikes rates and functional connectivity across '...
        'electrodes may be caused by the spikes themselves.</p>'],...
        pretty_exp_html(stats.medians(2)),pretty_exp_html(stats.medians(1)),...
        stats.Tpos,get_p_html(stats.pval),fig_name);
    
    annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
    annotation('textbox',[0.33 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
    annotation('textbox',[0.67 0.91 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
    
end

fclose(fid);
print(gcf,[plot_folder,fig_name],'-dpng')

%% Figure 5 - null model construction
if doing_from_github == 0 || doing_from_github == 2
    null_model_conceptual
end

%% Table 2 - model comparisons
model_comp_table

%% Figure 6 and Fig S5 - model
for ia = 1:2
    
    if ia == 1
        fname = 'results.html';
        model = brain_model;
        fig_name = 'Fig 6';
        fid = fopen([plot_folder,fname],'a');
        fprintf(fid,'<p><br><b>Predicting the SOZ</b></br>');
        fprintf(fid,['We next developed a classifier to predict whether an individual '...
        'electrode contact would be designated as belonging to the SOZ or not. To account for '...
        'expected spatial bias in electrode sampling (clinicians preferentially place '...
        'electrodes closer to the SOZ), we first constructed a null model, which provides '...
        'an estimate of the accuracy of predicting the SOZ based entirely on spatial location (Fig 5). '...
        'This is an estimate of how accurately we can localize the SOZ before even considering '...
        'EEG data.']);
    else
        fname = 'supplemental_results.html';
        model = aal_model;
        fig_name = 'Fig S5';
        fid = fopen([plot_folder,fname],'a');
        fprintf(fid,'<p><br><b>Predicting SOZ</b></br>');
        fprintf(fid,['We repeated the classification analysis using the AAL atlas '...
            'rather than the Brainnetome atlas.']);
    end
    
    
    


    model_info = model.all_out.model_info;
    mnames = {sprintf('Chance: AUC 0.50'),...
        sprintf('Null model: AUC %1.2f',mean(model_info(2).all_auc)),...
        sprintf('Connectivity + null: AUC %1.2f',mean(model_info(3).all_auc)),...
        sprintf('Spikes + null: AUC %1.2f',mean(model_info(4).all_auc)),...
        sprintf('All: AUC %1.2f',mean(model_info(5).all_auc))};
    mbars = {'-','--','-.','--o'};
    mleg = nan(3,1);
    figure
    set(gcf,'position',[380 147 833 650])
    tiledlayout(2,1,'tilespacing','tight','padding','tight');
    
    nexttile
    plot([0 1],[0 1],'k:','linewidth',2)
    hold on
    count = 0;
    for im = [2 3 4 5]
        count = count+1;
        if im == 2
            mname = 'Null model';
        elseif im == 3
            mname = 'Connectivity';
        elseif im == 4
            mname = 'Spikes';
        elseif im == 5
            mname = 'All';
        end
        if im == 5
            mleg(count) = line_fewer_markers(model_info(im).x,model_info(im).ym,10,mbars{count},'linewidth',2,...
                'color',def_colors(count,:));
        else
            mleg(count) = plot(model_info(im).x,model_info(im).ym,'linewidth',2,...
                'color',def_colors(count,:));
        end


    end
    set(gca,'fontsize',15)
    xlabel('False positive rate')
    ylabel('True positive rate')
    legend(mnames,'fontsize',15,'location','southeast')
    title('SOZ prediction')
    
    fprintf(fid,[' We compared the performance of the null model against a model including connectivity '...
        'data, a model including spike rate data, '...
        'and a model including both spikes and connectivity (%sA). This analysis included all patients '...
        'with electrode localizations and accurate spike detections (%d patients). The model '...
        'was trained on 2/3 of the patients, and tested on the remaining 1/3. 1,000 random '...
        'splits of patients into testing and training data were performed to obtain model statistics. '...
        'These results show that, first, '...
        'even a null model ignoring EEG data substantially outperforms a chance model (null model AUC = %1.2f). It also shows '...
        'incremental improvement in adding connectivity data (AUC = %1.2f; '...
        'albeit not as good as with adding spike data, AUC = %1.2f), with '...
        'still better performance when combining spike and connectivity data '...
        '(AUC = %1.2f; Table 2 shows statistics of model comparisons).</p>'],...
        fig_name,model.all_out.npts_remain,mean(model_info(2).all_auc),mean(model_info(3).all_auc),...
        mean(model_info(4).all_auc),mean(model_info(5).all_auc));
    
    % Loop over models
    nexttile
    count = 0;
    all_ps = nan(4,1);
    scale = 2;
    txticklocs = [];
    txticklabels = {};
    legp = nan(4,1);
    for im = [2 3 4 5]
        
        count = count + 1;
        % grab bootstrap aucs and stats
        aucs = model.all_out.stereo_vs_not.model(im).all;
        naucs = size(aucs,1);
        p = model.all_out.stereo_vs_not.model(im).p;
        all_ps(count) = p;
        
        % Loop over stereo vs not stereo
        for is = 1:2
            if is == 1
                sstereo = 'Stereo';
            else
                sstereo = 'Non-stereo';
            end
            % plot all aucs
            legp(count) = plot(im-1+(is-1.5)/scale+0.05*randn(naucs,1),aucs(:,is),'o',...
                'color',def_colors(count,:),'linewidth',2);
            hold on
            % plot mean auc
            plot([im-1+(is-1.5)/scale-0.1 im-1+(is-1.5)/scale+0.1],...
                [mean(aucs(:,is)) mean(aucs(:,is))],...
                'color',def_colors(count,:),'linewidth',2)
            txticklocs = [txticklocs im-1+(is-1.5)/scale];
            txticklabels = [txticklabels sstereo];
            
        end
        
        
    end
    
    % plot bars and p values
    yl = ylim;
    ybar = yl(1) + 1.05 * (yl(2)-yl(1));
    yp = yl(1) + 1.1 * (yl(2)-yl(1));
    ynew = [yl(1) yl(1) + 1.15 * (yl(2)-yl(1))];
    ylim(ynew)
    for im = 1:4
        plot([im-0.5/scale,im+0.5/scale],[ybar ybar],'k','linewidth',2)
        text(im,yp,get_p_text(all_ps(im)),'horizontalalignment','center','fontsize',15);
    end
    
    legend(legp,{'Null model','Connectivity','Spikes','All'},...
    'fontsize',15,'position',[0.5108 0.0730 0.1438 0.1192]);
    xticks(txticklocs)
    xticklabels(txticklabels)
    ylabel('Model AUC');
    set(gca,'fontsize',15)
    title('Model performance in stereo-EEG versus grid/strip/depth implantations')
    
    fprintf(fid,['<p>We anticipated that the spatial null model and other models ('...
        'all of which incorporate the spatial null model information) '...
        'might perform better in patients with stereo-EEG than in patients with '...
        'grid/strip/depth implantations because electrode spacing is uniform '...
        'in grids and strips. To test this, we compared the performance '...
        'of models trained and tested only on patients with stereo-EEG implantations '...
        '(%d patients who also had accurate spike detections and electrode localizations) against those '...
        'trained on patients with grid/strip/depth implantations (%d patients). As expected, the performance  '...
        'of each model was higher in the case of stereo-EEG implantations (%sB; Table S2).</p>'],...
        model.all_out.stereo_vs_not.n_stereo_remain,...
        model.all_out.stereo_vs_not.n_not_stereo_remain,fig_name);
    
    

    

    beta = model.all_out.glme_stuff.model(5).glm.Coefficients{5,2};
    se = model.all_out.glme_stuff.model(5).glm.Coefficients{5,3};
    t = model.all_out.glme_stuff.model(5).glm.Coefficients{5,4};
    p = model.all_out.glme_stuff.model(5).glm.Coefficients{5,6};
    or = exp(beta);
    ci95 = [exp(beta - 1.96*se),exp(beta + 1.96*se)];
    fprintf(fid,['<p>Finally, we used the estimate of the model coeffiencients to '...
        'assess how an electrodes''s connectivity is associated with the likelihood '...
        'of being in the SOZ, controlling for spatial bias and spike rates. For this analysis we used a logistic '...
        'mixed effects model trained on all patients. The patient identifier was a random effect. '...
        'Holding covariates constant, the odds of an electrode being a SOZ '...
        'electrode decreased by %1.1f%% (95%% CI OR [%1.2f-%1.2f]) for each additional '...
        'normalized connectivity unit (<i>t</i> = %1.1f, %s).'],(1-or)*100,ci95(1),...
        ci95(2),t,get_p_html(p));

    if ia == 1
        fprintf(fid,['The negative odds ratio implies that, controlling for spatial sampling '...
            'bias and spike rates, SOZ electrodes tend to have lower average connectivity.</p>']);
    end
    
    annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
    annotation('textbox',[0 0.40 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
    print(gcf,[plot_folder,fig_name],'-dpng')
    fclose(fid);
end

%% Table S2 - sEEG vs grid/strip/depth implantation
model_implant_table

close all

end