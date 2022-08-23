function subsample_elecs_test(brain_out,aal_out,plot_folder)



for ia = 1:2
    
    if ia == 1
        atlas_out = brain_out;
        fid = fopen([plot_folder,'supplemental_results.html'],'a');
        fprintf(fid,'<p><br><b>Testing effect of electrode number on symmetric coverage-Brainnetome atlas</b></br>');
        fig_name = 'Fig S2';
        
        fprintf(fid,['We next tested how the number of electrodes in each atlas region affected'...
            ' connectivity results in our symmetric coverage analysis, first using the Brainnetome atlas.']);
        
    else
        atlas_out = aal_out;
        fid = fopen([plot_folder,'supplemental_results.html'],'a');
        fprintf(fid,'<p><br><b>Testing effect of electrode number on symmetric coverage - AAL atlas</b></br>');
        fig_name = 'Fig S3';
        
        fprintf(fid,['We also tested how the number of electrodes in each atlas region affected'...
            ' connectivity results in our symmetric coverage analysis using the AAL atlas.']);
    end
    
    
    
    nb = length(atlas_out);

    % unpack brain_out - be careful not to confuse with aal results
    %unpack_any_struct(atlas_out);
    figure
    set(gcf,'position',[1 1 1400 900])
    tiledlayout(2,6,'tilespacing','tight','padding','tight')

    % Show bilateral coverage
    nexttile([1 3])
    
    %% Get single values for things that shouldn't change across iterations
    % THESE DO CHANGE!!! I THINK IT MAKES SENSE
    %{
    all_bilateral = atlas_out(1).all_bilateral;
    nsymmetric = atlas_out(1).n_coverage.nsymmetric;
    any_symmetric = atlas_out(1).n_coverage.any_symmetric;
    neither_lat = atlas_out(1).neither_lat;
    n_soz_all = atlas_out(1).n_analyses.n_soz_all;
    n_holo_hem = atlas_out(1).n_analyses.n_holo_hem;
    n_soz_intra = atlas_out(1).n_soz_intra;
    % These should be identical
    for ib = 1:nb
        assert(isequal(all_bilateral,atlas_out(ib).all_bilateral))
        assert(isequal(nsymmetric,atlas_out(ib).n_coverage.nsymmetric))
        assert(isequal(any_symmetric,atlas_out(ib).n_coverage.any_symmetric));
        assert(isequal(neither_lat,atlas_out(ib).neither_lat));
        assert(isequal(n_soz_all,atlas_out(ib).n_analyses.n_soz_all));
        assert(isequal(n_holo_hem,atlas_out(ib).n_analyses.n_holo_hem));
        assert(isequal(n_soz_intra,atlas_out(ib).n_soz_intra));
    end
    %}
    
    %%  Average other stuff across iterations
    all_all_bilateral = repmat(atlas_out(1).all_bilateral,1,1,nb);
    all_soz_non_soz_ordered_atlas = repmat(atlas_out(1).soz_non_soz_ordered_atlas,1,1,1,nb);
    all_soz_all = repmat(atlas_out(1).soz_all,1,1,nb);
    all_hemi = repmat(atlas_out(1).hemi,1,1,nb);
    all_soz_intra = repmat(atlas_out(1).soz_intra,1,1,nb);
    for ib = 1:nb
        all_soz_non_soz_ordered_atlas(:,:,:,ib) = atlas_out(ib).soz_non_soz_ordered_atlas;
        all_soz_all(:,:,ib) = atlas_out(ib).soz_all;
        all_hemi(:,:,ib) = atlas_out(ib).hemi;
        all_soz_intra(:,:,ib) = atlas_out(ib).soz_intra;
         all_all_bilateral(:,:,ib) = atlas_out(ib).all_bilateral;
    end
    soz_non_soz_ordered_atlas = nanmean(all_soz_non_soz_ordered_atlas,4);
    soz_all = nanmean(all_soz_all,3);
    hemi = nanmean(all_hemi,3);
    soz_intra = nanmean(all_soz_intra,3);
    all_bilateral = nanmean(all_all_bilateral,3);
    
    
    bilat_cov = all_bilateral(~atlas_out(1).neither_lat,:);
    left_bilat_cov = bilat_cov(atlas_out(1).left==1,:);
   % left_names = atlas_out.atlas_names(atlas_out.left == 1);
    pm = turn_nans_gray(left_bilat_cov);
    xticklabels([])
    yticks(1:size(left_bilat_cov,1))
    %yticklabels(left_names);
    yticklabels([])
    hold on
    xl = xlim;
    new_xl = [xl(1) - 0.15*(xl(2)-xl(1)),xl(2)];
    line_pos = xl(1) - 0.09*(xl(2)-xl(1));
    arrow_pos = xl(1) - 0.04*(xl(2)-xl(1));
    xlim(new_xl)
    colormap(gca,[0.5,0.5,0.5;1 1 1]);
    title('Regions with symmetric coverage')
    xlabel('Patient')
    ylabel('Region')
    set(gca,'fontsize',20)
    
    %% Also show some of the clusters
    if ia == 1
        % Show temporal neocortical cluster
        mid_pos = 47;
        text(line_pos,mid_pos,'Lateral\newlinetemporal','horizontalalignment','center',...
    'rotation',90,'fontsize',20)
    	text(arrow_pos,mid_pos,'\rightarrow','fontsize',20,'horizontalalignment','left','fontweight', 'bold')
        
        % Show mesial temporal cluster
        mid_pos = 108;
        text(line_pos,mid_pos,'Mesial\newlinetemporal','horizontalalignment','center',...
    'rotation',90,'fontsize',20)
        text(arrow_pos,mid_pos,'\rightarrow','fontsize',20,'horizontalalignment','left','fontweight', 'bold')
    else
        
        % Show mesial tempora;
        mid_pos = 19;
        text(line_pos,mid_pos,'Mesial\newlinetemporal','horizontalalignment','center',...
    'rotation',90,'fontsize',20)
        text(arrow_pos,mid_pos,'\rightarrow','fontsize',20,'horizontalalignment','left','fontweight', 'bold')
        
        % show lateral temporal
        mid_pos = 43;
        text(line_pos,mid_pos,'Lateral\newlinetemporal','horizontalalignment','center',...
    'rotation',90,'fontsize',20)
    	text(arrow_pos,mid_pos,'\rightarrow','fontsize',20,'horizontalalignment','left','fontweight', 'bold')
        
    end
        

  

    %% Show atlas
    
    
    
    nexttile([1 3])
    pretty_matrix(nanmean(soz_non_soz_ordered_atlas(~atlas_out(1).neither_lat,~atlas_out(1).neither_lat,:),3),...
        {'SOZ\newlinehemisphere','non-SOZ\newlinehemisphere'},sum(atlas_out(1).left),'Connectivity (r)',0);
    caxis(gca,[-1 1])
    title('Average connectivity (symmetric coverage only)')
    xlabel('Region')
    set(gca,'fontsize',20)
   
    % Show tests
    %need TO ADD N FOR EACH TEST
    nexttile([1 2])
    stats = paired_plot(soz_all,'Connectivity',{'to SOZ','to contralateral region'});
    title({'Connectivity','to SOZ vs. contralateral'})
    set(gca,'fontsize',20)
    fprintf(fid,['<p>We compared the average connectivity to the SOZ region with that to the '...
        'region contralateral to the SOZ. ']);

    fprintf(fid,[' The average connectivity to the SOZ region (median %s) was lower than that to '...
        ' the region contralateral to the SOZ (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sC).</p>'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name);

    nexttile([1 2])
    stats = paired_plot(hemi,'Intra-hemispheric connectivity',{'in SOZ side','in non-SOZ side'},1);
    title({'Hemispheric connectivity','ipsilateral vs. contralateral to SOZ'})
    set(gca,'fontsize',20)
    fprintf(fid,['<p>We next compared the intra-hemispheric connectivity between '...
        'the side of the SOZ and the contralateral hemisphere.']);

    fprintf(fid,[' The average intra-hemispheric connectivity on the side of the SOZ (median %s) was also lower than that in '...
        ' the contralateral hemisphere (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%sD).'...
        ' This suggests that the change in connectivity in epilepsy is broad and affects the entire hemisphere.</p>'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name);

    nexttile([1 2])
    stats = paired_plot(soz_intra,'Intrinsic connectivity',{'in SOZ','in contralateral region'},1);
    title({'Intrinsic connectivity','in SOZ vs contralateral region'})
    set(gca,'fontsize',20)
    fprintf(fid,['<p>We next compared the intrinsic connectivity within the SOZ regions '...
        'and that of the contralateral region.']);

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


    fclose(fid);
end


end
