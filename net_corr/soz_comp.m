function soz_comp(all)

soz = all.net.ccep.is_soz;
ccep = all.net.ccep.data;

% reduce soz
new_labels = all.net.ccep.orig_labels;
lia = ismember(new_labels,all.labels);
new_labels = new_labels(lia);
soz = soz(lia);

% re-arrange
[~,locb] = ismember(all.labels,new_labels);
new_labels = new_labels(locb);
soz = soz(locb);

if ~isequal(new_labels,all.labels), error('why'); end

% Get outdegree and indegree
outdegree = (nansum(ccep,1))';
indegree = nansum(ccep,2);

% Compare for soz and non soz
figure
tiledlayout(1,2)

nexttile
plot(1+0.05*rand(sum(soz==1),1),outdegree(soz),'o')
hold on
plot(2+0.05*rand(sum(soz==0),1),outdegree(~soz),'o')
xlim([0 3])
xticks([1 2])
xticklabels({'SOZ','non-SOZ'});
ylabel('Outdegree')

nexttile
plot(1+0.05*rand(sum(soz==1),1),indegree(soz),'o')
hold on
plot(2+0.05*rand(sum(soz==0),1),indegree(~soz),'o')
xlim([0 3])
xticks([1 2])
xticklabels({'SOZ','non-SOZ'});
ylabel('Indegree')

end