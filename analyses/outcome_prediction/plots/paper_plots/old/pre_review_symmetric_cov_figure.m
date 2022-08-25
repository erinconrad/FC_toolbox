function symmetric_cov_figure(brain_out,aal_out,plot_folder,corr_out)

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
        fig_name = 'Fig S1';
        
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
        {'Left','Right'},sum(atlas_out.left),[],1);
    colormap(gca,[0.5,0.5,0.5;1 1 1]);
    title('Regions with symmetric coverage')
    xlabel('Patient')
    ylabel('Region')
    set(gca,'fontsize',20)

    fprintf(fid,[' Across patients, there were on average %1.1f (range %d-%d) anatomical regions '...
        ' (%1.1f when separately counting left and right) with bilateral electrode coverage. '...
        '%d of %d patients had regions with bilateral electrode coverage (%sA).'],...
        mean(atlas_out.n_coverage.nsymmetric),min(atlas_out.n_coverage.nsymmetric),max(atlas_out.n_coverage.nsymmetric),...
        mean(atlas_out.n_coverage.nsymmetric)*2,atlas_out.n_coverage.any_symmetric,sum(corr_out.pts_with_any_locs),fig_name);

    % Show atlas
    nexttile([1 3])
    pretty_matrix(nanmean(atlas_out.soz_non_soz_ordered_atlas(~atlas_out.neither_lat,~atlas_out.neither_lat,:),3),...
        {'SOZ\newlinehemisphere','non-SOZ\newlinehemisphere'},sum(atlas_out.left),'Connectivity (r)',0);
    caxis(gca,[-1 1])
    title('Average connectivity (symmetric coverage only)')
    xlabel('Region')
    set(gca,'fontsize',20)
    fprintf(fid,[' Visually, in a network averaged across patients, network edges in the hemisphere opposite the SOZ were '...
        'stronger than those in the hemisphere of the SOZ (%sB, brighter colors indicate higher connectivity).</p>'],fig_name);

    % Show tests
    %need TO ADD N FOR EACH TEST
    nexttile([1 2])
    stats = paired_plot(atlas_out.soz_all,'Connectivity',{'to SOZ','to contralateral region'});
    title({'Connectivity','to SOZ vs. contralateral'})
    set(gca,'fontsize',20)
    fprintf(fid,['<p>We compared the average connectivity to the SOZ with that to the '...
        'region contralateral to the SOZ. We included patients who had electrode coverage of both the SOZ '...
        'and the contralateral region in this analysis (%d patients).'],...
        atlas_out.n_analyses.n_soz_all);

    fprintf(fid,[' The average connectivity to the SOZ (median %s) was lower than that to '...
        ' the region contralateral to the SOZ (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sC).</p>'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name);

    nexttile([1 2])
    stats = paired_plot(atlas_out.hemi,'Intra-hemispheric connectivity',{'in SOZ side','in non-SOZ side'},1);
    title({'Hemispheric connectivity','ipsilateral vs. contralateral to SOZ'})
    set(gca,'fontsize',20)
    fprintf(fid,['<p>We next compared the intra-hemispheric connectivity between '...
        'the side of the SOZ and the contralateral hemisphere. For this analysis we included '...
        'patients with unilateral epilepsy and at least two symmetrically-implanted atlas regions in '...
        'each hemisphere, which allowed us to calculate intra-hemispheric connectivity (%d patients).'],...
        atlas_out.n_analyses.n_holo_hem);

    fprintf(fid,[' The average intra-hemispheric connectivity on the side of the SOZ (median %s) was also lower than that in '...
        ' the contralateral hemisphere (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sD).'...
        ' This suggests that the change in connectivity in epilepsy is broad and affects the entire hemisphere.</p>'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name);

    nexttile([1 2])
    stats = paired_plot(atlas_out.soz_intra,'Intrinsic connectivity',{'in SOZ','in contralateral region'},1);
    title({'Intrinsic connectivity','in SOZ vs contralateral region'})
    set(gca,'fontsize',20)
    fprintf(fid,['<p>We next compared the intrinsic connectivity within the SOZ '...
        'and that of the contralateral region. For this analysis we included '...
        'patients whose SOZ spanned at least two atlas regions, and who had electrode coverage '...
        'of these regions and of the contralateral regions (%d patients).'],...
        atlas_out.n_soz_intra);

    fprintf(fid,[' The intrinsic connectivity within the SOZ '...
        '(between one SOZ region and other SOZ regions) (median %s) was no different from the intrinsic connectivity '...
        ' within the same regions in the contralateral hemisphere (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sE).'...
        ' This suggests no focal within-SOZ change in connectivity within the limits of this analysis.'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name);

    % Add annotations
    annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
    annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',30,'linestyle','none')
    annotation('textbox',[0 0.48 0.1 0.1],'String','C','fontsize',30,'linestyle','none')
    annotation('textbox',[0.33 0.48 0.1 0.1],'String','D','fontsize',30,'linestyle','none')
    annotation('textbox',[0.68 0.48 0.1 0.1],'String','E','fontsize',30,'linestyle','none')

    print(gcf,[plot_folder,fig_name],'-dpng')
    
    if ia == 1
        fprintf(fid,[' Results were similar when using the AAL rather than the Brainnetome '...
            'atlas for parcellating brain regions (Supplemental Results; Fig S1).']);
    end
    
    
    fprintf(fid,[' Results when studying coherence rather than Pearson correlation networks '...
        'were more heterogeneous, seen for different frequency bands in the two different atlases '...
        '(Supplemental Results; Fig S2).</p>']);
    fclose(fid);
end


end
