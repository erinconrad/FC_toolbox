function out = outcome_plot_orders(things,sozs)

myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];
grayColor = 0.75*[1 1 1];
markersize = 3;
npts = length(things);


%[all_ranks,all_soz_ranks,nchance,all,successes] = get_ranks(things,sozs,rates,which,min_rate);
[all_ranks,all_soz_ranks,~,soz_chance] = simple_rate_rank(things,sozs);

% Alternate statistical test
%{
In this test, I see whether the median rate is higher for SOZ than chance
for each patient (whereas in the prior I see if the median rate RANKING is
higher for SOZ than chance). I think this is better because when I take
ranks matlab will force equal things to have different ranks, and so this
seems more accurate. The p values are slightly different but both <0.001
%}
if 1
    assert(isequal(any(isnan(soz_chance),2),cellfun(@(x) sum(x)==0,sozs)))
    
    % remove ties (these shouldn't count in significance testing, right?)
    % and nans
    nans = any(isnan(soz_chance),2);
    ties = soz_chance(:,1) == soz_chance(:,2);
    nans_or_ties = nans | ties;
    soz_chance(nans_or_ties,:) = [];
    
    npts_alt = size(soz_chance,1);
 
    alt_successes = soz_chance(:,1) > soz_chance(:,2);
    alt_failures = soz_chance(:,1) < soz_chance(:,2);
    
    
    
    assert(sum(alt_successes) + sum(alt_failures) == npts_alt);
    
    if sum(alt_successes==0) < sum(alt_successes==1)
        pval_binom_alt = 2*binocdf(sum(alt_successes==0),npts_alt,0.5);
    else
        pval_binom_alt = 2*binocdf(sum(alt_successes==1),npts_alt,0.5);
    end
    
end

%% Double check with another binomial test % remove when publish code as requires someone else's code
% alt alt
pout=myBinomTest(sum(alt_successes==0),npts_alt,0.5,'two');

% check it gives the same result
assert(abs(pout-pval_binom_alt)<0.01);

%% Re-order by # of electrodes
num_elecs = cellfun(@length,all_ranks);

% remove those with no elecs
no_elecs = num_elecs == 0;
all_ranks(no_elecs) = [];
all_soz_ranks(no_elecs) = [];
%successes(no_elecs) = [];

num_elecs = cellfun(@length,all_ranks);
npts = length(num_elecs);

[num_elecs,I] = sort(num_elecs);
all_ranks = all_ranks(I);
all_soz_ranks = all_soz_ranks(I);


%% Plot stuff
for i = 1:npts
    plot(i,all_ranks{i},'o','color',grayColor,'markersize',markersize);
    hold on
    sp = plot(i,nanmedian(all_soz_ranks{i}),'p',...
        'markerfacecolor',[0.2660 0.6040 0.2880],'markeredgecolor',[0.2660 0.6040 0.2880],...
        'linewidth',2,'markersize',markersize+5);
    cp = plot(i,median(1:num_elecs(i)),'o','markeredgecolor',[0.3 0.3 0.3],...
        'markerfacecolor',[0.25 0.25 0.25],...
        'linewidth',2,'markersize',markersize+2);
    
end



% Get median across all patients
all_all_soz_ranks = [];
for i = 1:npts
    all_all_soz_ranks = [all_all_soz_ranks;(all_soz_ranks{i})'];
    
end
median_rank = median(all_all_soz_ranks);

%out.n = length(successes);
%out.nsuc = sum(successes==1);
out.pval = pval_binom_alt;
out.median_rank = median_rank;

xlim([0 npts])
xl = xlim;
xbar = xl(1) + 1.02*(xl(2)-xl(1));
xtext = xl(1) + 1.04*(xl(2)-xl(1));
newxl = [xl(1) xl(1) + 1.05*(xl(2)-xl(1))];
plot([xbar xbar],[1 num_elecs(end)],'k-','linewidth',2)
if pval_binom_alt >= 0.05
text(xtext-1,(1+num_elecs(end))/2,get_asterisks(pval_binom_alt,1),'rotation',90,...
        'horizontalalignment','center','fontsize',16)

else
    text(xtext,(1+num_elecs(end))/2,get_asterisks(pval_binom_alt,1),'rotation',90,...
        'horizontalalignment','center','fontsize',20)
end
xlim(newxl)
legend([sp,cp],{'Median SOZ', 'Chance'},'fontsize',15,'location','northwest')

end