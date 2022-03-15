function neurology_figure5



%% Parameters
nb = 1000;
plot_type = 'scatter';
nblocks = 6;
myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];



locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
%out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/neurology/'];

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');
fprintf(fid,'<p><b>Localizing the seizure onset zone with spikes</b><br>');

figure
set(gcf,'position',[10 10 900 600])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

%% Get stuff
rate_sw = out.bin_out.all_elecs_rates_sw;
rate = out.bin_out.all_elecs_rates;
soz = out.bin_out.all_is_soz;
npts = length(soz);

%% SOZ spike rate ranking
nexttile([1 2])
stats_out = plot_orders(rate,soz);
hold on
xticklabels([])
xlabel('Patient')
ylabel('Electrode spike rate rank')
set(gca,'fontsize',15)
title('Seizure onset zone - spike rate ranking')

fprintf(fid,['To determine whether sleep disproportionately increased spikes '...
    'in the SOZ relative to non-SOZ electrodes, we examined the ranking of electrodes'...
    ' by spike rates. We first examined whether SOZ electrodes '...
    'had a higher spike rate than expected by chance. The median SOZ electrode '...
    'rank across patients was %1.1f, implying that across patients, %d '...
    'non-SOZ electrodes tended to have higher spike rates than the SOZ electrodes. '...
    'The median SOZ rank was higher (closer to 1) than the median overall electrode rank '...
    'in %d of %d patients, which is more than expected by chance '...
    '(Binomial test, %s) (Fig. 5A). These findings imply that the SOZ has '...
    'more frequent spikes than expected by chance.'],stats_out.median_rank,stats_out.median_rank-1,...
    stats_out.nsuc,...
    stats_out.n,get_p_html(stats_out.pval));


%% Spike rate ranking sleep vs wake

% For sleep/wake ranking comparison, generate set in which any nans in
% sleep or wake are removed (to compare same set of electrodes when
% comparing ranks in wake and asleep)
any_nans = cellfun(@(x) any(isnan(x),2),rate_sw,'uniformoutput',false);
npts = length(any_nans);
rate_sw_nan_removed = rate_sw;
soz_nan_removed = soz;
rate_nan_removed = rate;
for i = 1:npts
    curr_nans = any_nans{i};
    curr_rate_sw = rate_sw{i};
    curr_soz = soz{i};
    
    % remove nans
    rate_sw_nan_removed{i} = curr_rate_sw(~curr_nans,:);
    soz_nan_removed{i} = curr_soz(~curr_nans);
    rate_nan_removed{i} = curr_soz(~curr_nans);
end

% confirm no nans left for wake or sleep
assert(sum(cellfun(@(x) any(isnan(x),'all'),rate_sw_nan_removed)) == 0);

% Get wake and sleep rate
wake_rate = cellfun(@(x) x(:,1), rate_sw_nan_removed,'uniformoutput',false);
sleep_rate = cellfun(@(x) x(:,2), rate_sw_nan_removed,'uniformoutput',false);


[wake_all_ranks,wake_soz_ranks] = simple_rate_rank(wake_rate,soz_nan_removed);
wake_soz_ranks = cellfun(@nanmedian,wake_soz_ranks);

[sleep_all_ranks,sleep_soz_ranks] = simple_rate_rank(sleep_rate,soz_nan_removed);
sleep_soz_ranks = cellfun(@nanmedian,sleep_soz_ranks);

soz_rank_sw_rate = [wake_soz_ranks,sleep_soz_ranks];

% Plot it
nexttile([1 1])
stats = plot_paired_data(soz_rank_sw_rate',{'wake','sleep'},sprintf('Spike rate rank'),'paired','scatter','ranking');
title(sprintf('SOZ spike rate ranking'))

% Results text
fprintf(fid,[' The spike rate ranking of SOZ electrodes was higher '...
    ' (closer to 1) in sleep (median rank %1.1f)'...
    ' than in wake (median rank %1.1f) '...
    '(Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 5B).'...
    ' This implies that sleep preferentially increases spike rates '...
    'in the SOZ relative to other electrodes.</p>'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));

% confirm that the total number of electrodes being ranked is the same for
% wake and sleep (otherwise unfair comparison)
assert(isequal(cellfun(@length,wake_all_ranks),cellfun(@length,sleep_all_ranks)))

%% Fancy model
% Do it once for the purpose of generating my ORs
soz_roc_out = classify_soz(0);
%roc = soz_roc_out.roc;
%auc = soz_roc_out.auc;
sleep_or = soz_roc_out.sleep_or;
wake_or = soz_roc_out.wake_or;
sleep_ci95 = soz_roc_out.sleep_ci95;
wake_ci95 = soz_roc_out.wake_ci95;
sleep_t= soz_roc_out.sleep_t;
wake_t= soz_roc_out.wake_t;
sleep_p= soz_roc_out.sleep_p;
wake_p= soz_roc_out.wake_p;

sleep_or = (sleep_or - 1)*100;
wake_or = (wake_or - 1)*100;


all_roc = nan(nb,1000,2);
all_auc = nan(nb,1);
all_disc = nan(nb,1);
all_soz_class = cell(nb,1);
all_non_soz_class = cell(nb,1);
% Now do it 1,000 times!!!
for ib = 1:nb
    soz_roc_out = classify_soz(1);
    all_roc(ib,:,:) = soz_roc_out.roc;
    all_auc(ib) = soz_roc_out.auc;
    all_disc(ib) = soz_roc_out.disc;
    all_soz_class{ib} = soz_roc_out.class_soz;
    all_non_soz_class{ib} = soz_roc_out.class_no_soz;
end

nexttile([1 1])
median_auc = nanmedian(all_auc);
auc_ci = [prctile(all_auc,2.5) prctile(all_auc,97.5)];

% to get the shaded bits, need to do fancy interpolation to get the x's
% into the same space
[unify_x,unify_y] = unify_roc(all_roc);

ym = median(unify_y,1);


ly = prctile(unify_y,25,1);
uy = prctile(unify_y,75,1);


%% ROC text
fprintf(fid,['<p>Finally, we tested how accurately spike rates could '...
    'classify electrodes as SOZ versus non-SOZ. We trained a logistic regression'...
    ' classifier on a randomly chosen two-thirds of the patients, using the normalized '...
    'electrode-specific average spike rate in each of the sleep and wake states as predictor variables '...
    'and the binary classification of the electrode as SOZ or non-SOZ as '...
    'the response variable. Holding spike rate in wake constant, the odds'...
    ' of an electrode being a SOZ electrode increased by %1.1f%% (95%% CI [%1.2f-'...
    '%1.2f]) for each additional normalized spike rate in sleep unit '...
    '(<i>t</i> = %1.1f, %s). Holding spike rate in sleep constant,'...
    ' the odds of an electrode being a SOZ electrode increased by %1.1f%% (95%% CI [%1.2f-'...
    '%1.2f]) for each additional normalized spike rate in wake unit ('...
    '<i>t</i> = %1.1f, %s). This larger odds ratio associated with spikes in sleep provides additional support that spikes in sleep localize the SOZ better than spikes in wake.'...
    ' Applying this model to the remaining'...
    ' one-third of patients held as testing data accurately '...
    'classified electrodes as SOZ versus non-SOZ '...
    '(Median AUC across 1,000 random splits of testing and training data = %1.2f (95%% CI [%1.2f-%1.2f]), Fig. 5C).'...
    ' These results imply that a simple model incorporating just spike rates and sleep/wake state can accurately localize the SOZ.</p>'],...
    sleep_or,sleep_ci95(1),sleep_ci95(2),sleep_t,get_p_html(sleep_p),wake_or,...
    wake_ci95(1),wake_ci95(2),wake_t,get_p_html(wake_p),median_auc,auc_ci(1),auc_ci(2));

%% Build confusion matrix
disc = median(all_disc);
all_mat = nan(2,2,nb);
all_sens = nan(nb,1);
all_spec = nan(nb,1);
all_ppv = nan(nb,1);
all_npv = nan(nb,1);
for ib = 1:nb
    tout = sleep_conf_matrix(all_soz_class{ib},all_non_soz_class{ib},disc);
    all_mat(:,:,ib) = tout.mat;
    all_sens(ib) = tout.sensitivity;
    all_spec(ib) = tout.specificity;
    all_ppv(ib) = tout.ppv;
    all_npv(ib) = tout.npv;
end

%% Alt confusion matrix to optimize PPV
disc_opt = 0.2;
new_all_mat = nan(2,2,nb);
new_all_sens = nan(nb,1);
new_all_spec = nan(nb,1);
new_all_ppv = nan(nb,1);
new_all_npv = nan(nb,1);
for ib = 1:nb
    tout = sleep_conf_matrix(all_soz_class{ib},all_non_soz_class{ib},disc_opt);
    new_all_mat(:,:,ib) = tout.mat;
    new_all_sens(ib) = tout.sensitivity;
    new_all_spec(ib) = tout.specificity;
    new_all_ppv(ib) = tout.ppv;
    new_all_npv(ib) = tout.npv;
end

%% Add text
%{
fprintf(fid,[' As an example of how this classifier would perform at predicting the SOZ,'...
    ' the point on the median ROC curve that best maximizes the sum of sensitivity and specificity'...
    ' results in a median sensitivity of %1.1f%% and specifity of %1.1f%% for detecting the SOZ across the 1,000 random test datasets.'...
    ' This corresponds to a positive predictive value (the probability that a predicted SOZ is a true SOZ) of %1.1f%%'...
    ' and a negative predictive value of %1.1f%%. An alternate value on the ROC curve optimizing specificity at the'...
    ' expense of sensitivity results in a positive predictive value of %1.1f%% and a negative predictive value of %1.1f%%.'...
    ' These results imply that a simple model incorporating just spike rates and sleep/wake state can accurately localize the SOZ.</p>'],...
    median(all_sens)*100,median(all_spec)*100,median(all_ppv)*100,median(all_npv)*100,...
    median(new_all_ppv)*100,median(new_all_npv)*100);
%}

mp = shaded_error_bars_fc(unify_x,ym,[ly',uy'],'k');
%plot(roc(:,1),roc(:,2),'k-','linewidth',2)
hold on
plot([0 1],[0 1],'k--','linewidth',2)

% Add the example points
point1 = [1-median(all_spec),median(all_sens)];
point2 = [1-median(new_all_spec),median(new_all_sens)];

% snap the example points to the graph
point1 = find_closest_point(point1,[unify_x',ym']);
point2 = find_closest_point(point2,[unify_x',ym']);

%{
p1=plot(point1(1),point1(2),'d','markeredgecolor',[0.4940, 0.1840, 0.5560],...
    'markerfacecolor',[0.4940, 0.1840, 0.5560],'markersize',15);
p2 = plot(point2(1),point2(2),'o','markeredgecolor',[0.3010, 0.7450, 0.9330],...
    'markerfacecolor',[0.3010, 0.7450, 0.9330],'markersize',15);
%}
legend(sprintf('Median AUC %1.2f',median_auc),'location','southeast','fontsize',15,'box','off')
xlabel('False positive rate')
ylabel('True positive rate')
%{
legend([mp,p1,p2],...
    {sprintf('Median AUC %1.2f',median_auc),sprintf('PPV = %1.1f%%, NPV = %1.1f%%',...
    median(all_ppv)*100,median(all_npv)*100),sprintf('PPV = %1.1f%%, NPV = %1.1f%%',...
    median(new_all_ppv)*100,median(new_all_npv)*100)},...
    'location','southeast','fontsize',15,'box','off')
%}
set(gca,'fontsize',15)
title('SOZ identification accuracy')

%% Add annotations
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.415 0.1 0.1],'String','B','fontsize',20,'linestyle','none')
annotation('textbox',[0.5 0.415 0.1 0.1],'String','C','fontsize',20,'linestyle','none')


fclose(fid);
print([out_folder,'fig5'],'-depsc')
close(gcf)


end

function [unify_x,unify_y] = unify_roc(all_roc)

nx = 10000;
unify_x = linspace(0,1,nx);
unify_y = nan(size(all_roc,1),nx);

for i = 1:size(all_roc,1)
    curr_x = all_roc(i,:,1);
    for ix = 1:nx
        curr_ix = unify_x(ix);
        
        % find the curr x closest to this ix
        [~,I] = min(abs(curr_ix-curr_x));
        
        % fill the y in with the corresponding y
        corr_y = all_roc(i,I,2);
        unify_y(i,ix) = corr_y;
        
    end
end

end


function closest = find_closest_point(test,points)

dist = vecnorm(points-test,2,2);
[~,closest] = min(dist);
closest = points(closest,:);

end
