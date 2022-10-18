function epilepsia_supplemental_range_threshold(scores,soz,X,ym)


range_thresh = [.03 0.1 0.2 0.4];
nthresh = length(range_thresh);

figure
set(gcf,'position',[10 10 900 1000])
tiledlayout(4,2,'tilespacing','tight','padding','tight')
for t = 1:nthresh
    desired_threshold = range_thresh(t);

    [npv,ppv,~,~,mat,~,acc,sens,spec] = cellfun(@(x,y,z) individual_threshold_stats(x,y,desired_threshold),...
    scores,soz,'uniformoutput',false);
    sens = cell2mat(sens);
    spec = cell2mat(spec);
    npv = cell2mat(npv);
    ppv = cell2mat(ppv);
    acc = cell2mat(acc);
    mat = cat(3,mat{:}); mat = mat./nansum(mat,[1 2]);
    all_mat = nanmean(mat,3);
    
    % combine mats
    %{
    all_mat = zeros(2,2);
    for ip = 1:length(mat)
        if any(isnan(mat{ip}),'all')
            continue
        else
            all_mat = all_mat + mat{ip};
        end
    end
    %}
    % combine mats

    
    % Calculate true positives, etc.
    tp = all_mat(2,2); % pred true and actual true
    tn = all_mat(1,1); % predicted false and actual false
    fn = all_mat(2,1); % predicted false and actual true
    fp = all_mat(1,2); % predicted true and actual false
    
    alt_total_ppv = tp/(fp+tp);
    alt_total_npv = tn/(tn+fn);
    alt_total_acc = (tp+tn)/(tp+tn+fp+fn);
    all_sens = tp/(tp+fn);
    all_spec = tn/(tn+fp);

    alt_alt_sens = mat(2,2,:)./(mat(2,2,:)+mat(2,1,:));

    % plot roc
    nexttile
    plot(X,ym,'linewidth',2)
    hold on 
    plot([0 1],[0 1],'k--','linewidth',2)
    plot(1-nanmean(spec),nanmean(sens),'*','markersize',15,'linewidth',2);
    set(gca,'fontsize',15)
    xlabel('False positive rate')
    ylabel('True positive rate')
    title('SOZ classification accuracy')

    nexttile
    turn_nans_gray([1 0;0 1])
    colormap(gca,[0.8000, 0.420, 0.42;0.5, 0.75, 0.5])
    xticks(1:2)
    xticklabels({'Not SOZ','SOZ'})
    yticks(1:2)
    yticklabels({'Not SOZ','SOZ'})
    xlabel('Predicted')
    ylabel('Actual')
    hold on
    conf_text = {'True negatives','False positives';'False negatives','True positives'};
    for ic = 1:2
        for jc = 1:2
            text(ic,jc,sprintf('%s:\n%1.1f%%',conf_text{jc,ic},all_mat(jc,ic)*100),'horizontalalignment','center','fontsize',20,'fontweight','bold')
        end
    end
    title(sprintf(['Threshold %1.2f,'...
         ' PPV: %1.1f%%, NPV: %1.1f%%'],...
    desired_threshold,...
    nanmean(ppv)*100,nanmean(npv)*100))
    set(gca,'fontsize',15)

end





end