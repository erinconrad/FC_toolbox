function model_figure(brain_model,aal_model,plot_folder,doing_from_github)

def_colors = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.9290, 0.6940, 0.1250;...
    0.4940, 0.1840, 0.5560];

for ia = 1:2
    
    if ia == 1
        fname = 'results.html';
        model = brain_model;
        fig_name = 'Fig 5';
        fid = fopen([plot_folder,fname],'a');
        fprintf(fid,'<p><br><b>Predicting the SOZ</b></br>');
        fprintf(fid,['We next developed a classifier to predict whether an individual '...
        'electrode contact belongs to the SOZ or not. To account for '...
        'expected spatial bias in electrode sampling (clinicians preferentially place '...
        'electrodes closer to the SOZ), we first constructed a null model, which provides '...
        'an estimate of the accuracy of predicting the SOZ electrodes based entirely on spatial location (Fig 4). '...
        'This is an estimate of how accurately we can localize the SOZ electrodes before even considering '...
        'EEG data.']);
        tab_name = '2';
    else
        fname = 'supplemental_results.html';
        model = aal_model;
        fig_name = 'Fig S6';
        fid = fopen([plot_folder,fname],'a');
        fprintf(fid,'<p><br><b>Predicting SOZ - AAL atlas</b></br>');
        fprintf(fid,['We repeated the classification analysis using the AAL atlas '...
            'rather than the Brainnetome atlas.']);
        tab_name = 'S1';
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
            if doing_from_github == 1
                mleg(count) = plot(model_info(im).x,model_info(im).ym,mbars{count},'linewidth',2,...
                    'color',def_colors(count,:));
            else
                mleg(count) = line_fewer_markers(model_info(im).x,model_info(im).ym,10,mbars{count},'linewidth',2,...
                    'color',def_colors(count,:));
            end
        else
            mleg(count) = plot(model_info(im).x,model_info(im).ym,'linewidth',2,...
                'color',def_colors(count,:));
        end


    end
    set(gca,'fontsize',20)
    xlabel('False positive rate')
    ylabel('True positive rate')
    legend(mnames,'fontsize',20,'location','southeast')
    title('SOZ prediction')
    
    prctile95null = prctile(model_info(2).all_auc,[2.5 97.5]);
    fprintf(fid,[' We compared the performance of the null model against a model including connectivity '...
        'data, a model including spike rate data, '...
        'and a model including both spikes and connectivity (%sA). This analysis included all patients '...
        'with electrode localizations and accurate spike detections (%d patients). The model '...
        'was trained on 2/3 of the patients, and tested on the remaining 1/3. 1,000 random '...
        'splits of patients into testing and training data were performed to obtain model statistics. '...
        'These results show that, first, '...
        'even a null model based on electrode placement, ignoring EEG data substantially outperforms a chance model (null model AUC = %1.2f, 95%% CI %1.2f-%1.2f). It also shows '...
        'incremental improvement in adding connectivity data (AUC = %1.2f; '...
        'albeit not as good as with adding spike data, AUC = %1.2f), with '...
        'still better performance when combining spike and connectivity data '...
        '(AUC = %1.2f; Table %s shows statistics of model comparisons). Notably, using a conservative '...
        'two-sided bootstrap method to compare AUCs of different models, a model '...
        'adding functional connectivity data did not significantly outperform the spatial '...
        'null model alone (%s).</p>'],...
        fig_name,model.all_out.npts_remain,...
        mean(model_info(2).all_auc),prctile95null(1),prctile95null(2),...
        mean(model_info(3).all_auc),...
        mean(model_info(4).all_auc),mean(model_info(5).all_auc),tab_name,get_p_html(model.all_out.model_p(3,2)));
    
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
    'fontsize',20,'position',[0.5108 0.120 0.1438 0.1192]);
    legend boxoff
    xticks(txticklocs)
    xticklabels(txticklabels)
    ylabel('Model AUC');
    set(gca,'fontsize',20)
    title('Model performance in stereo-EEG versus grid/strip/depth implantations')
    
    fprintf(fid,['<p>We anticipated that the spatial null model and other models ('...
        'all of which incorporate the spatial null model information) '...
        'might perform better in patients with stereo-EEG than in patients with '...
        'grid/strip/depth implantations because electrode spacing is uniform '...
        'in grids and strips. To test this, we compared the performance '...
        'of models trained and tested only on patients with stereo-EEG implantations '...
        '(%d patients who also had accurate spike detections and electrode localizations) against those '...
        'trained on patients with grid/strip/depth implantations (%d patients). The performance  '...
        'of each model was higher in the case of stereo-EEG implantations, although the distributions heavily overlapped,'...
        ' and no difference was significant (%sB; Table S2).</p>'],...
        model.all_out.stereo_vs_not.n_stereo_remain,...
        model.all_out.stereo_vs_not.n_not_stereo_remain,fig_name);
    
    

    

    beta_mean = model.all_out.glme_stuff.stats.mean;
    ci95_beta = model.all_out.glme_stuff.stats.CI_95;
    p = model.all_out.glme_stuff.stats.p;
    or = exp(beta_mean);
    ci95 = [exp(ci95_beta(1)),exp(ci95_beta(2))];
    fprintf(fid,['<p>Finally, we used the estimate of the model coefficients to '...
        'assess how an electrode''s connectivity is associated with the likelihood '...
        'of being in the SOZ, controlling for spatial bias and spike rates. For this analysis we used a logistic '...
        'mixed effects model trained on all patients. The patient identifier was a random effect. '...
        'Holding covariates constant, the odds of an electrode being a SOZ '...
        'electrode decreased by %1.1f%% (bootstrap 95%% CI OR [%1.2f-%1.2f]) for each additional '...
        'normalized connectivity unit (bootstrap %s).'],(((1-or)*100)),(ci95(1)),...
        (ci95(2)),get_p_html(p));

    if ia == 1
        fprintf(fid,[' The odds ratio less than 1 implies that, controlling for spatial sampling '...
            'bias and spike rates, a lower average connectivity increases the likelihood of an electrode being in the SOZ.</p>']);
    end
    
    annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
    annotation('textbox',[0 0.40 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
    print(gcf,[plot_folder,fig_name],'-dpng')
    fclose(fid);
end


end