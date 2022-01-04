function out = plot_orders(rates,sozs)

myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];
grayColor = 0.75*[1 1 1];
markersize = 2;
npts = length(rates);


%[all_ranks,all_soz_ranks,nchance,all,successes] = get_ranks(things,sozs,rates,which,min_rate);
[all_ranks,all_soz_ranks,successes] = simple_rate_rank(rates,sozs);

%% Re-order by # of electrodes
num_elecs = cellfun(@length,all_ranks);
[num_elecs,I] = sort(num_elecs);
all_ranks = all_ranks(I);
all_soz_ranks = all_soz_ranks(I);


%% Plot stuff
for i = 1:npts
    plot(i,all_ranks{i},'o','color',grayColor,'markersize',markersize);
    hold on
    sp = plot(i,nanmedian(all_soz_ranks{i}),'*','color',myColours(1,:),'linewidth',2,'markersize',markersize+2);
    cp = plot(i,median(1:num_elecs(i)),'*','color',myColours(2,:),'linewidth',2,'markersize',markersize+2);
    
end


pval_binom = 2*binocdf(sum(successes==0),length(successes),0.5);

% Get median across all patients
all_all_soz_ranks = [];
for i = 1:npts
    all_all_soz_ranks = [all_all_soz_ranks;(all_soz_ranks{i})'];
    
end
median_rank = median(all_all_soz_ranks);

out.n = length(successes);
out.nsuc = sum(successes==1);
out.pval = pval_binom;
out.median_rank = median_rank;

xlim([0 npts])
xl = xlim;
xbar = xl(1) + 1.03*(xl(2)-xl(1));
xtext = xl(1) + 1.06*(xl(2)-xl(1));
newxl = [xl(1) xl(1) + 1.08*(xl(2)-xl(1))];
plot([xbar xbar],[1 num_elecs(end)],'k-','linewidth',2)
if pval_binom >= 0.05
text(xtext-1,(1+num_elecs(end))/2,get_asterisks(pval_binom,1),'rotation',90,...
        'horizontalalignment','center','fontsize',16)

else
    text(xtext,(1+num_elecs(end))/2,get_asterisks(pval_binom,1),'rotation',90,...
        'horizontalalignment','center','fontsize',20)
end
xlim(newxl)
legend([sp,cp],{'Median SOZ', 'Chance'},'fontsize',15,'location','northwest')

end