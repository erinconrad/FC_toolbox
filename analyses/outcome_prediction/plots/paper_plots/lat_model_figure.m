function lat_model_figure(lat_info,plot_folder)
which_models = 1:3;

main_fig_name = 'Fig 3';
supp_fig_name = 'Fig S4';
    

for ia = 1:2
    
    figure
    set(gcf,'position',[380 147 1400 450])
    tiledlayout(1,length(which_models),'tilespacing','tight','padding','tight');
    
    if ia == 1
        fid = fopen([plot_folder,'results.html'],'a');
        atlas_text = 'Brainnetome';
        fig_name = main_fig_name;
        fprintf(fid,'<p><br><b>Epilepsy lateralization</b></br>');
        fprintf(fid,['We asked how well intra-hemispheric functional connectivity '...
            'could lateralize epilepsy.']);
    else
        atlas_text = 'AAL';
        fid = fopen([plot_folder,'supplemental_results.html'],'a');
        fig_name = supp_fig_name;
        fprintf(fid,'<p><br><b>Epilepsy lateralization - AAL atlas</b></br>');
        fprintf(fid,'We repeated the laterality model using the AAL atlas rather than the Brainnetome atlas.');
    end
    unpack_any_struct(lat_info(ia))

    for im = which_models
        nexttile
    
        if im == 1
            mname = 'Connectivity';
        elseif im == 2
            mname = 'Spikes';
        elseif im == 3
            mname = 'Connectivity + Spikes';
        end
        
        turn_nans_gray([1 0;0 1])
        colormap(gca,[0.8500, 0.3250, 0.0980;0, 0.4470, 0.7410])
        xticks(1:2)
        xticklabels({'Left','Right'})
        yticks(1:2)
        yticklabels({'Left','Right'})
        xlabel('Predicted')
        ylabel('True')
        hold on
        for ic = 1:2
            for jc = 1:2
                text(ic,jc,sprintf('%d',all_conf_table(jc,ic,im)),'horizontalalignment','center','fontsize',25,'fontweight','bold')
            end
        end
        title(sprintf('%s %s\nAccuracy: %1.1f%%\nPPV: %1.1f%%, NPV: %1.1f%%',atlas_text,mname,...
            all_accuracy(im)*100,...
            all_lpv(im)*100,all_rpv(im)*100))
        set(gca,'fontsize',20)

    end
    
    %% Annotations
    annotation('textbox',[0 0.90 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
    annotation('textbox',[0.35 0.90 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
    annotation('textbox',[0.68 0.90 0.1 0.1],'String','C','fontsize',25,'linestyle','none')

  
    print(gcf,[plot_folder,fig_name],'-dpng');
    
    fprintf(fid,[' We restricted analysis to patients with unilateral SOZ '...
    'and at least two symmetrically-implanted atlas regions in '...
    'each hemisphere, which allowed us to calculate intra-hemispheric connectivity (%d patients).'],...
    lat_info(ia).removed.noriginal-lat_info(ia).removed.nremoved_non_unilateral-lat_info(ia).removed.nremoved_missing_fc);

    fprintf(fid,[' The average accuracy of SOZ laterality predictions in the test patients was %1.1f%%'...
        ' (positive predictive value (PPV) for detecting left-sided seizure onset: %1.1f%%, negative predictive value (NPV): %1.1f%%) (%sA).'],...
        lat_info(ia).all_accuracy(1)*100,lat_info(ia).all_lpv(1)*100,lat_info(ia).all_rpv(1)*100,fig_name);
    
    if ia == 1
        fprintf(fid,[' This model '...
        'performed similarly to one using spike rates rather than connectivity (%sB). '...
        'Using a model based on the AAL atlas rather than the Brainnetome atlas '...
        'yielded similar results (%s).</p>'],fig_name,supp_fig_name);
    else
        fprintf(fid,[' This model '...
        'had somewhat better performance than one using spike rates rather than connectivity (%sB).</p>'],supp_fig_name);
    end
    
    
    
end






close all
fclose(fid);

end