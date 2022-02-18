function plot_orders_mats(things,sozs)

myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];
grayColor = 0.75*[1 1 1];
markersize = 5;
npts = size(things,2);
nregions = size(things,2);



%% Get median thing for SOZ and for all non-nans
soz_chance = nan(npts,2);
for i = 1:npts
    % curr things and soz
    curr_things = things(:,i);
    curr_soz = sozs(:,i);
    
    % remove nans from things and soz
    any_nan = isnan(curr_things);
    curr_things(any_nan) = [];
    curr_soz(any_nan) = [];
    
    % get median for each
    soz_chance(i,:) = [median(curr_things(curr_soz)) median(curr_things)];
    
end

% remove nan rows
soz_chance(any(isnan(soz_chance),2),:) = [];
npts = size(soz_chance,1);
successes = soz_chance(:,1) < soz_chance(:,2);
pval_binom = 2*binocdf(sum(successes==0),sum(~isnan(successes)),0.5);

[nranks,median_soz_ranks,alt_successes] = simple_rate_rank_mat(things,sozs,successes);


%% Re-order by # of electrodes
num_elecs = nranks;

% remove those with no elecs
no_elecs = num_elecs == 0;
median_soz_ranks(no_elecs) = [];
nranks(no_elecs) = [];

num_elecs = nranks;
npts = length(num_elecs);

[num_elecs,I] = sort(num_elecs);
median_soz_ranks = median_soz_ranks(I);
nranks = nranks(I);


%% Plot stuff
for i = 1:npts
    if nranks(i) == 0, continue; end
    plot(i,1:nranks(i),'o','color',grayColor,'markersize',markersize);
    hold on
    sp = plot(i,median_soz_ranks(i),'*','color',myColours(1,:),'linewidth',2,'markersize',markersize+2);
    cp = plot(i,median(1:nranks(i)),'*','color',myColours(2,:),'linewidth',2,'markersize',markersize+2);
    
end


%pval_binom = 2*binocdf(sum(successes==0),sum(~isnan(successes)),0.5);

out.n = length(successes);
out.nsuc = sum(successes==1);
out.pval = pval_binom;

xlim([0 npts])
xl = xlim;
xbar = xl(1) + 1.02*(xl(2)-xl(1));
xtext = xl(1) + 1.04*(xl(2)-xl(1));
newxl = [xl(1) xl(1) + 1.05*(xl(2)-xl(1))];
plot([xbar xbar],[1 nranks(end)],'k-','linewidth',2)
if pval_binom >= 0.05
text(xtext-1,(1+nranks(end))/2,get_asterisks(pval_binom,1),'rotation',90,...
        'horizontalalignment','center','fontsize',16)

else
    text(xtext,(1+nranks(end))/2,get_asterisks(pval_binom,1),'rotation',90,...
        'horizontalalignment','center','fontsize',20)
end
xlim(newxl)
legend([sp,cp],{'Median SOZ', 'Chance'},'fontsize',15,'location','northwest')

end