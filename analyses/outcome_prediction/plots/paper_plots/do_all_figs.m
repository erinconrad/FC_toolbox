function do_all_figs

close all

%% Parameters
freqs = {'Delta (0.5-4 Hz)','Theta (4-8 Hz)','Alpha (8-12 Hz)','Beta (12-30 Hz)','Gamma (30-80 Hz)'};
def_colors = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.9290, 0.6940, 0.1250;...
    0.4940, 0.1840, 0.5560];

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

%% Delete the current results file
if exist([plot_folder,'results.html'])~=0
    delete([plot_folder,'results.html']);
end

if exist([plot_folder,'supplemental_results.html'])~=0
    delete([plot_folder,'supplemental_results.html']);
end

%% Do all analyses
% Load results from model
brain_model = load([model_folder,'model_stuff_brainnetome.mat']);
aal_model = load([model_folder,'model_stuff_aal_bernabei.mat']);

% symmetric coverage tests for both atlases
aal_out = load([model_folder,'symm_cov_aal_bernabei.mat']);
aal_out = aal_out.nout;
brain_out = load([model_folder,'symm_cov_brainnetome.mat']);
brain_out = brain_out.nout;

% Spike-fc corr
corr_out = spike_fc_correlation;

% density out
dens_out = load([model_folder,'dens_model.mat']);
dens_out = dens_out.out;

%% Fig 1 - conceptual fig
main_conceptual_figure

%% Supplemental fig 1 - density model construction
ip = 80;
all_fc = dens_out.all_fc;
all_locs = dens_out.all_locs;
vec_dens = dens_out.vec_dens;
vec_conn = dens_out.vec_conn;
resid = dens_out.resid;
sr = dens_out.sr;
g = dens_out.g;

    
figure
set(gcf,'position',[10 10 1000 700])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

% Functional connectivity
nexttile
conn = all_fc{ip};
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
D = make_interdist_matrix(all_locs{ip});
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

nexttile
%[xbins,ybins] = bin_data_get_means(vec_dens,vec_conn,1e3);
plot(vec_dens,vec_conn,'o')
hold on
temp_x = [min(vec_dens):0.1*(max(vec_dens)-min(vec_dens)):max(vec_dens)];
%temp_y = (f.p1 * temp_x + f.p2)./(temp_x + f.q1);
temp_y = (g.p1 + temp_x*g.p2);
plot(temp_x,temp_y,'k','linewidth',3)
ylim([-1 1])
xlim([min(vec_dens) max(vec_dens)])
xlabel('Density (mm)')
ylabel('Functional connectivity (r)')
set(gca,'fontsize',15)
title({'Density-functional connectivity correlation','(all patients)'})

nexttile
conn = resid{ip};
turn_nans_gray(conn)
c = colorbar;
%caxis([-1 1])
ylabel(c,'Normalized correlation (model residuals)')
set(gca,'fontsize',15)
xticklabels([])
yticklabels([])
xlabel('Electrode')
ylabel('Electrode')
title({'Density-normalized correlation','(single patient)'})

annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.42 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.42 0.1 0.1],'String','D','fontsize',25,'linestyle','none')

print(gcf,[plot_folder,'FigS1'],'-dpng')

%% Figure 2 - symmetric coverage test for brainnetome and Supplemental Fig 2 (same but AAL)
% NEED TO THINK ABOUT N
% Prep section in text
for ia = 1:2
    
    if ia == 1
        atlas_out = brain_out;
        fid = fopen([plot_folder,'results.html'],'a');
        fprintf(fid,'<p><br><b>Connectivity with symmetric coverage constraint</b></br>');
        fig_name = 'Fig2';
        
        fprintf(fid,['To control for spatial sampling bias, we first compared connectivity '...
        'in the SOZ against that in the contralateral region while restricting analysis '...
        'to regions with symmetric electrode coverage.']);
        
    else
        atlas_out = aal_out;
        fid = fopen([plot_folder,'supplemental_results.html'],'a');
        fprintf(fid,'<p><br><b>Connectivity with symmetric coverage constraint - AAL atlas</b></br>');
        fig_name = 'FigS2';
        
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
        '%d of %d patients had any regions with bilateral electrode coverage (Fig. 2A).'],...
        mean(atlas_out.n_coverage.nsymmetric),min(atlas_out.n_coverage.nsymmetric),max(atlas_out.n_coverage.nsymmetric),...
        mean(atlas_out.n_coverage.nsymmetric)*2,atlas_out.n_coverage.any_symmetric,atlas_out.total_n_for_symmetric);
    % why is n different?

    % Show atlas
    nexttile([1 3])
    pretty_matrix(nanmean(atlas_out.soz_non_soz_ordered_atlas(~atlas_out.neither_lat,~atlas_out.neither_lat,:),3),...
        {'SOZ\newlinehemisphere','non-SOZ\newlinehemisphere'},sum(atlas_out.left),'r',0);
    caxis(gca,[-1 1])
    title('Average connectivity (symmetric coverage only)')
    xlabel('Region')

    fprintf(fid,[' Visually, in a network averaged across patients, network edges in the hemisphere opposite the SOZ were '...
        'stronger than those in the hemisphere of the SOZ (Fig. 2B, brighter colors indicate higher connectivity).']);

    % Show tests
    %need TO ADD N FOR EACH TEST
    nexttile([1 2])
    stats = paired_plot(atlas_out.soz_all,'Average connectivity',{'SOZ','contralateral region','contralateral region'});
    title('Connectivity to SOZ vs. contralateral region')
    
    fprintf(fid,[' We compared the average connectivity to the SOZ with that to the '...
        'region contralateral to the SOZ. We included patients who had electrode coverage of both the SOZ '...
        'and the contralateral region in this analysis (%d patients).'],...
        atlas_out.n_analyses.n_soz_all);

    fprintf(fid,[' The average connectivity to the SOZ (median %s) was lower than that to '...
        ' the region contralateral to the SOZ (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 2C).'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval));

    nexttile([1 2])
    stats = paired_plot(atlas_out.hemi,'Intra-hemispheric connectivity',{'SOZ side','non-SOZ side','non-SOZ side'});
    title({'Hemispheric connectivity','ipsilateral vs. contralateral to SOZ'})
    
    fprintf(fid,[' We next compared the intra-hemispheric connectivity between '...
        'the side of the SOZ and the contralateral hemisphere. For this analysis we included '...
        'patients with unilateral epilepsy and at least two symmetrically-implanted atlas regions in '...
        'each hemisphere, which allowed us to calculate intra-hemispheric connectivity (%d patients).'],...
        atlas_out.n_analyses.n_holo_hem);

    fprintf(fid,[' The average intra-hemispheric connectivity on the side of the SOZ (median %s) was also lower than that in '...
        ' the contralateral hemisphere (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 2D).'...
        ' This suggests that the change in connectivity in epilepsy is broad and affects the entire hemisphere.'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval));

    nexttile([1 2])
    stats = paired_plot(atlas_out.soz_intra,'Intrinsic connectivity',{'SOZ','contralateral region','contralateral region'});
    title({'Intrinsic connectivity','in SOZ vs contralateral region'})
    
    fprintf(fid,[' We next compared the intrinsic connectivity within the SOZ '...
        'and that of the contralateral region. For this analysis we included '...
        'patients whose SOZ spanned at least two atlas regions, and who had electrode coverage '...
        'of these regions and of the contralateral regions (%d patients).'],...
        atlas_out.n_soz_intra);

    fprintf(fid,[' The intrinisic connectivity within the SOZ '...
        '(between one SOZ region and other SOZ regions) (median %s) was no different from the intrinsic connectivity '...
        ' within the same regions in the contralateral hemisphere (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 2E).'...
        ' This suggests no focal within-SOZ change in connectivity within the limits of this analysis.'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval));

    % Add annotations
    annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
    annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',30,'linestyle','none')
    annotation('textbox',[0 0.44 0.1 0.1],'String','C','fontsize',30,'linestyle','none')
    annotation('textbox',[0.32 0.44 0.1 0.1],'String','D','fontsize',30,'linestyle','none')
    annotation('textbox',[0.66 0.44 0.1 0.1],'String','E','fontsize',30,'linestyle','none')

    print(gcf,[plot_folder,fig_name],'-dpng')
    
    if ia == 1
        fprintf(fid,[' Results were similar when using the AAL rather than the Brainnetome '...
            'atlas for parcellating brain regions (Supplemental Results, Fig. S2).</p>']);
    end
    
    fprintf(fid,[' Results for studying coherence rather than Pearson correlation networks '...
        'were more heterogeneous, seen for different frequency bands in the two different atlases '...
        '(Supplemental Results, Fig. S3).</p>']);
    fclose(fid);
end

%% Supplemental Fig 3 - coherence (both atlases)

figure
set(gcf,'position',[1 1 1440 900])
tiledlayout(2,5,'tilespacing','tight','padding','tight')
fig_name = 'FigS3';
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
        stats = paired_plot(squeeze(soz_coh_all(:,f,:)),'Coherence',{'SOZ','contralateral region','contralateral region'});
        title(freqs{f})
        if stats.pval < 0.05
            fprintf(fid,[' The average %s coherence to the SOZ (median %s) was lower than that to '...
                ' the region contralateral to the SOZ (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 1C).'],...
                freqs{f},pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval));
        else
            fprintf(fid,[' There was no significant difference between the %s coherence to the SOZ (median %s) and that to '...
                ' the region contralateral to the SOZ (median %1.1e) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 1C).'],...
                freqs{f},pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval));
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
        fig_name = 'Fig3';
        fid = fopen([plot_folder,'results.html'],'a');
        fprintf(fid,'<p><br><b>Epilepsy lateralization</b></br>');
        fprintf(fid,['We asked how well intra-hemispheric functional connectivity '...
            'could lateralize epilepsy.']);
    else
        atlas_out = aal_out;  
        fig_name = 'FigS4';
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
                    'at lateralizing the SOZ (%s) (Fig. 3).'],conf_out.accuracy*100,...
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
                    '(Supplemental Results, Fig. S4).']);
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
set(gcf,'position',[1 100 370 1000])
tiledlayout(1,3,'tilespacing','compact','padding','tight');
fid = fopen([plot_folder,'results.html'],'a');
fprintf(fid,'<p><br><b>Correlation between functional connectivity and spikes</b></br>');

fprintf(fid,['We next measured the correlation between functional connecitivity and spike rates. '...
    'For each patient, we took the Spearman rank correlation between the average spike rate '...
    'for each electrode and the average node strength for each electrode.']);
for i = 1
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
    mean_corr = nanmean(thing);
    xlabel('Patient')
    xticklabels([])
    ylabel('Spike rate - connectivity correlation (\rho)')
    title('Spike-functional connectivity correspondence')
    [~,p,~,stats] = ttest(thing);
    legend(pp,sprintf('mean \\rho = %1.2f, %s',mean_corr,get_p_text(p)),'location','southeast','fontsize',15)
    set(gca,'fontsize',15)
    
    nexttile
    paired_plot(all_all_chs_corr,'Average connectivity',{'before spike','during spike'});
    
    nexttile
    paired_plot(all_sp_chs_corr,'Average connectivity',{'before spike','during spike'});
    
end

fprintf(fid,[' The mean correlation was &rho; = %1.2f, which was significantly less than zero '...
    '(<i>t</i>-test, <i>t</i>(%d) = %1.2f, %s (Fig. 4).'],mean_corr,stats.df,stats.tstat,get_p_html(p));
fprintf(fid,[' This implies that electrodes with higher spike rates tend to be less connected '...
    'to other electrodes.</p>']);
fclose(fid);

%annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
%annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
print(gcf,[plot_folder,'Fig4'],'-dpng')

%% Figure 5 - null model construction
null_model_conceptual

%% Table 2 - model comparisons
model_comp_table

%% Figure 6 - model
for ia = 1:2
    
    if ia == 1
        fname = 'results.html';
        model = brain_model;
        fig_name = 'Fig6';
        fid = fopen([plot_folder,fname],'a');
        fprintf(fid,'<p><br><b>Predicting SOZ</b></br>');
        fprintf(fid,['We next developed a classifier to predict whether an individual '...
        'electrode contact would be designated as belonging to the SOZ or not. Given that electrode '...
        'sampling is spatially biased, we first constructed a null model, which provides '...
        'an estimate of the accuracy of predicting the SOZ based entirely on spatial location (Fig. 5). '...
        'This is an estimate of how accurately we can localize the SOZ before even considering '...
        'EEG data.']);
    else
        fname = 'supplemental_results.html';
        model = aal_model;
        fig_name = 'FigS5';
        fid = fopen([plot_folder,fname],'a');
        fprintf(fid,'<p><br><b>Predicting SOZ</b></br>');
        fprintf(fid,['We repeated the classification analysis using the AAL atlas '...
            'rather than the Brainnetome atlas (Fig. 5).']);
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
    print(gcf,[plot_folder,fig_name],'-dpng')

    fprintf(fid,[' We compared the performance of this null model against a model including connectivity '...
        'data, a model including spike rate data, '...
        'and a model including both spikes and connectivity (Fig. 6). The model '...
        'was trained on 2/3 of the patients, and tested on the remaining 1/3. 1,000 random '...
        'splits of patients into testing and training data were performed to obtain model statistics. '...
        'These results show that, first, '...
        'even a null model ignoring EEG data substantially outperforms a chance model. It also shows '...
        'incremental improvement in adding connectivity data (albeit not as good as with adding spike data), with '...
        'still better performance when using both spikes and connectivity data '...
        '(Table 2 shows statistics of model comparisons). </p>']);

    beta = model.all_out.glme_stuff.model(5).glm.Coefficients{5,2};
    se = model.all_out.glme_stuff.model(5).glm.Coefficients{5,3};
    t = model.all_out.glme_stuff.model(5).glm.Coefficients{5,4};
    p = model.all_out.glme_stuff.model(5).glm.Coefficients{5,6};
    or = exp(beta);
    ci95 = [exp(beta - 1.96*se),exp(beta + 1.96*se)];
    fprintf(fid,['<p>Finally, we used the estimate of the model coeffiencients to '...
        'assess how an electrodes''s connectivity is associated with the likelihood '...
        'of being in the SOZ, controlling for spatial bias and spike rates. For this analysis we used a generalized linear '...
        'mixed effects model trained on all patients. The patient identifier was a random effect. '...
        'Holding covariates constant, the odds of an electrode being a SOZ '...
        'electrode decreased by %1.1f%% (95%% CI OR [%1.2f-%1.2f]) for each additional '...
        'normalized node strength unit (<i>t</i> = %1.1f, %s).'],(1-or)*100,ci95(1),...
        ci95(2),t,get_p_html(p));

    if ia == 1
        fprintf(fid,['The negative odds ratio implies that, controlling for spatial sampling '...
            'bias and spike rates, SOZ electrodes tend to have lower average connectivity.</p>']);
    end

    fclose(fid);
end

close all

end