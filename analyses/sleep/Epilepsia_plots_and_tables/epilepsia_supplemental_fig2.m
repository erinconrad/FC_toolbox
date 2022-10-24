function epilepsia_supplemental_fig2

%% Parameters
min_rate = 0.1;
plot_type = 'scatter';
nblocks = 6;
%{
myColours = [0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250];
%}

myColours = [0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];



locations = fc_toolbox_locs;
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/epilepsia/'];

%% Prep output text file
fid = fopen([out_folder,'supplemental_results.html'],'a');
fprintf(fid,'<p><b>Alternate division of localization</b>');

%% Figure
figure
set(gcf,'position',[10 10 900 350])
tiledlayout(1,2,'tilespacing','tight','padding','tight')

%%  localization
loc = circ_out.all_locs;
mt = contains(loc,'mesial temporal');
nc = contains(loc,'cortical') | contains(loc,'other cortex'); % temporal neocortical and other cortex



%% Sleep-wake spike difference
rate_sw = bin_out.all_rates;
rate_diff = (rate_sw(:,2)-rate_sw(:,1));

if 0
   table(loc,mt,nc) 
end

[p,~,stats] = ranksum(rate_diff(mt),rate_diff(nc));
W = stats.ranksum;
nt = sum(~isnan(rate_diff(mt)));
ne = sum(~isnan(rate_diff(nc)));


% Calculate effect size
z = stats.zval;
n = nt+ne;
r = abs(z)/sqrt(n);


nexttile
plot(1+randn(sum(mt),1)*0.05,rate_diff(mt),'o','linewidth',2,'color',myColours(1,:))
hold on
plot([0.7 1.3],[nanmedian(rate_diff(mt)) nanmedian(rate_diff(mt))],...
    'linewidth',2,'color',myColours(1,:))
plot(2+randn(sum(nc),1)*0.05,rate_diff(nc),'o','linewidth',2,'color',myColours(2,:))
plot([1.7 2.3],[nanmedian(rate_diff(nc)) nanmedian(rate_diff(nc))],...
    'linewidth',2,'color',myColours(2,:))
xticks([1 2])
xticklabels({'Mesial temporal','Neocortical'})
ylabel('Sleep-wake spikes/elec/min')
title('Sleep-wake spike rate difference')
set(gca,'fontsize',15);
xlim([0 3])
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,sprintf('%s, effect size r = %1.2f',get_p_text(p),r),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)


ns = min([sum(mt),sum(nc)]);
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);

fprintf(fid,[' There was no significant difference in sleep-wake spike rate difference between patients'...
    ' with mesial temporal lobe epilepsy (median = %1.1f spikes/elecs/min) and patients with neocortical'...
    ' epilepsy (median = %1.1f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>mesial temporal</sub></i> = %d, <i>N<sub>neocortical</sub></i> = %d) ='...
    ' %1.1f, %s, effect size r = %1.2f) (Fig. 2H).'],nanmedian(rate_diff(mt)),nanmedian(rate_diff(nc)),...
    nt,ne,U,get_p_html(p),r);



%% Postictal rate change
nexttile([1 1])
all_pts_spikes_bins = sz_out.all_pts_spikes_bins;
nbins = size(all_pts_spikes_bins,2);
pre = 1:nbins/2;
post = nbins/2+1:nbins;
pre_post = [nanmean(all_pts_spikes_bins(:,pre),2),nanmean(all_pts_spikes_bins(:,post),2)];
rate_diff = (pre_post(:,2)-pre_post(:,1));
[p,~,stats] = ranksum(rate_diff(mt),rate_diff(nc));
W = stats.ranksum;
nt = sum(~(isnan(rate_diff(mt))));
ne = sum(~(isnan(rate_diff(nc))));
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;

% Calculate effect size
z = stats.zval;
n = nt+ne;
r = abs(z)/sqrt(n);

plot(1+randn(sum(mt),1)*0.05,rate_diff(mt),'o','linewidth',2,'color',myColours(1,:))
hold on
plot([0.7 1.3],[nanmedian(rate_diff(mt)) nanmedian(rate_diff(mt))],...
    'linewidth',2,'color',myColours(1,:))
plot(2+randn(sum(nc),1)*0.05,rate_diff(nc),'o','linewidth',2,'color',myColours(2,:))
plot([1.7 2.3],[nanmedian(rate_diff(nc)) nanmedian(rate_diff(nc))],...
    'linewidth',2,'color',myColours(2,:))
xticks([1 2])
xticklabels({'Mesial temporal','Neocortical'})
ylabel('Post-pre spikes/elec/min')
title('Post-preictal spike rate difference')
set(gca,'fontsize',15);
xlim([0 3])
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,sprintf('%s, effect size r = %1.2f',get_p_text(p),r),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)

% text


U = min([U1,U2]);
fprintf(fid,[' Patients with mesial temporal epilepsy also had a greater '...
    'pre-to-postictal increase in spike rates '...
    '(median = %1.2f spikes/elecs/min) than those with neocortical'...
    ' epilepsy (median = %1.2f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>mesial temporal</sub></i> = %d, <i>N<sub>neocortical</sub></i> = %d) ='...
    ' %1.1f, %s, effect size r = %1.2f).</p>'],nanmedian(rate_diff(mt)),nanmedian(rate_diff(nc)),...
    nt,ne,U,get_p_html(p),r);



%% Add annotations
annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.52 0.9 0.1 0.1],'String','B','fontsize',25,'linestyle','none')




print([out_folder,'FigS2'],'-depsc')
close(gcf)


end

function prc = prc_asleep(x)

prc = 100 * sum(x==0)/(sum(x==1)+sum(x==0));

end