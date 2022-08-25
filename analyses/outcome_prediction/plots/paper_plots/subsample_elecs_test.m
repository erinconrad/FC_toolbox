function subsample_elecs_test(brain_out,aal_out,plot_folder)

fig_name = 'Fig S2';
figure
set(gcf,'position',[1 1 1400 900])
tiledlayout(2,3,'tilespacing','compact','padding','tight')
subnames = {'A','B','C','D','E','F'};
count =0;

for ia = 1:2
    
    if ia == 1
        atlas_out = brain_out;
        fid = fopen([plot_folder,'supplemental_results.html'],'a');
        fprintf(fid,'<p><br><b>Testing effect of electrode number on symmetric coverage-Brainnetome atlas</b></br>');
        tit_text = 'Brainnetome';
        
        fprintf(fid,['We next tested how the number of electrodes in each atlas region affected'...
            ' connectivity results in our symmetric coverage analysis, first using the Brainnetome atlas.']);
        
    else
        atlas_out = aal_out;
        fid = fopen([plot_folder,'supplemental_results.html'],'a');
        fprintf(fid,'<p><br><b>Testing effect of electrode number on symmetric coverage - AAL atlas</b></br>');        
        fprintf(fid,['We also tested how the number of electrodes in each atlas region affected'...
            ' connectivity results in our symmetric coverage analysis using the AAL atlas.']);
        tit_text = 'AAL';
    end
    
    
    
    nb = length(atlas_out);

    
    
    %%  Average other stuff across iterations
    all_soz_all = repmat(atlas_out(1).soz_all,1,1,nb);
    all_hemi = repmat(atlas_out(1).hemi,1,1,nb);
    all_soz_intra = repmat(atlas_out(1).soz_intra,1,1,nb);
    for ib = 1:nb
        all_soz_all(:,:,ib) = atlas_out(ib).soz_all;
        all_hemi(:,:,ib) = atlas_out(ib).hemi;
        all_soz_intra(:,:,ib) = atlas_out(ib).soz_intra;
    end
    soz_all = nanmean(all_soz_all,3);
    hemi = nanmean(all_hemi,3);
    soz_intra = nanmean(all_soz_intra,3);
    
    fprintf(fid,[' We subsampled electrodes contributing to each region for the purpose '...
        'of connectivity measurements so that the number of electrodes was the same as the '...
        'number in the contralateral region. We performed this random subsampling '...
        '100 times and averaged the final connectivity measurements prior to performing '...
        'Wilcoxon signed rank tests comparing the SOZ to the contralateral region.']);
   
    % Show tests
    nexttile
    count = count+1;
    stats = paired_plot(soz_all,'Connectivity',{'to SOZ','\nto contralateral region'});
    title({'Connectivity','to SOZ vs. contralateral'},{[tit_text,' atlas']})
    set(gca,'fontsize',20)
    fprintf(fid,[' We compared the average connectivity to the SOZ region with that to the '...
        'region contralateral to the SOZ. ']);

    fprintf(fid,[' The average connectivity to the SOZ region (median %s) was lower than that to '...
        ' the region contralateral to the SOZ (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%s%s).'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name,subnames{count});

    nexttile
    count = count+1;
    stats = paired_plot(hemi,'Intra-hemispheric connectivity',{'in SOZ side','in non-SOZ side'},1);
    title({'Hemispheric connectivity','ipsilateral vs. contralateral to SOZ'},{[tit_text,' atlas']})
    set(gca,'fontsize',20)
    fprintf(fid,[' We next compared the intra-hemispheric connectivity between '...
        'the side of the SOZ and the contralateral hemisphere.']);

    fprintf(fid,[' The average intra-hemispheric connectivity on the side of the SOZ (median %s) was also lower than that in '...
        ' the contralateral hemisphere (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%s%s).'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name,subnames{count});

    nexttile
    count = count+1;
    stats = paired_plot(soz_intra,'Intrinsic connectivity',{'in SOZ','in contralateral region'},1);
    title({'Intrinsic connectivity','in SOZ vs contralateral region'},{[tit_text,' atlas']})
    set(gca,'fontsize',20)
    fprintf(fid,[' We next compared the intrinsic connectivity within the SOZ regions '...
        'and that of the contralateral region.']);

    fprintf(fid,[' The intrinsic connectivity within the SOZ '...
        '(between one SOZ region and other SOZ regions) (median %s) was no different from the intrinsic connectivity '...
        ' within the same regions in the contralateral hemisphere (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%s%s).'],...
        pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name,subnames{count});

    if ia==1
        fprintf(fid,'</p>');
    end
    
end

fprintf(fid,[' Overall, the similarity of these results to the primary analysis in which we did not subsample electrodes within a region '...
    'suggests that the number of electrodes contributing to a given '...
    'region does not have a large influence on the average inter-regional connectivity measurements.<p>']);

% Add annotations
annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
annotation('textbox',[0.32 0.9 0.1 0.1],'String','B','fontsize',30,'linestyle','none')
annotation('textbox',[0.68 0.9 0.1 0.1],'String','C','fontsize',30,'linestyle','none')
annotation('textbox',[0 0.4 0.1 0.1],'String','D','fontsize',30,'linestyle','none')
annotation('textbox',[0.32 0.4 0.1 0.1],'String','E','fontsize',30,'linestyle','none')
annotation('textbox',[0.68 0.4 0.1 0.1],'String','F','fontsize',30,'linestyle','none')

print(gcf,[plot_folder,fig_name],'-dpng')


fclose(fid);


end
