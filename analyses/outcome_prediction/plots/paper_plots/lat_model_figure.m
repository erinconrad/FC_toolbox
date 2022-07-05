function lat_model_figure(lat_info,plot_folder)

main_fid = fopen([plot_folder,'results.html'],'a');
figure
set(gcf,'position',[380 147 1000 450])
tiledlayout(1,2,'tilespacing','tight','padding','tight');
    

for ia = 1:2
    
    if ia == 1
        atlas_text = 'Brainnetome';
    else
        atlas_text = 'AAL';
    end
    unpack_any_struct(lat_info(ia))

    mnames = {sprintf('Connectivity'),...
        sprintf('Spikes'),...
        sprintf('Connectivity')};
    
    for im = 1
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
        title(sprintf('%s\nAccuracy: %1.1f%%\nPPV: %1.1f%%, NPV: %1.1f%%',atlas_text,...
            all_accuracy(im)*100,...
            all_lpv(im)*100,all_rpv(im)*100))
        set(gca,'fontsize',20)

    end
      
    
    print(gcf,[plot_folder,'Fig4'],'-dpng');
   
    
        
    

        
    
end

fid = main_fid;
fprintf(fid,'<p><br><b>Epilepsy lateralization</b></br>');
fprintf(fid,['We asked how well intra-hemispheric functional connectivity '...
    'could lateralize epilepsy.']);


fprintf(fid,[' We restricted analysis to patients with unilateral epilepsy '...
    'and at least two symmetrically-implanted atlas regions in '...
    'each hemisphere, which allowed us to calculate intra-hemispheric connectivity (%d patients).'],...
    removed.noriginal-removed.nremoved_non_unilateral-removed.nremoved_missing_fc);

fprintf(fid,[' We constructed a logistic regression classifier to predict seizure onset '...
    'lateralization. We used leave-one-out '...
    'cross-validation, separately keeping each patient as testing data while '...
    'training the model on the remaining patients. Predictors were '...
    'average left intra-hemispheric connectivity and average right intra-hemispheric connectivity.'...
    ' The average accuracy of seizure onset zone laterality predictions was %1.1f%% (Fig. 4A).'...
    ' Using a model based on the AAL atlas rather than the Brainnetome atlas '...
    'yielded an average accuracy of %1.1f%% (Fig. 4B).</p>'],...
    lat_info(1).all_accuracy(1)*100,...
    lat_info(2).all_accuracy(1)*100);

close all
fclose(fid);

end