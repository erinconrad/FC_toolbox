function epilepsia_supplemental_range_threshold(scores,soz,X,ym)

%% Locations
locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/epilepsia/'];

range_thresh = [.03 0.1 0.2 0.4];
nthresh = length(range_thresh);

myColours = [0 33 87;...
122 80 113;...    
227 124 29;...
    86 152 163]/255;

figure
set(gcf,'position',[10 10 1000 1000])
tiledlayout(4,2,'tilespacing','tight','padding','compact')
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
    plot(X,ym,'linewidth',2,'color',myColours(1,:))
    hold on 
    plot([0 1],[0 1],'k--','linewidth',2)
    plot(1-nanmean(spec),nanmean(sens),'*','markersize',15,'linewidth',2,'color',[163 2 52]/255);
    set(gca,'fontsize',20)
    xlabel('FPR')
    ylabel('TPR')
    title('SOZ classification accuracy')

    nexttile
    turn_nans_gray([1 0;0 1])
    %colormap(gca,[0.8000, 0.420, 0.42;0.5, 0.75, 0.5])
    colormap(gca,[206 128 128;161 197 203]/255)
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
    set(gca,'fontsize',20)

end
fontname(gcf,"calibri");

print(gcf,[out_folder,'FigS4'],'-depsc')
print(gcf,[out_folder,'FigS4'],'-dpng')



end