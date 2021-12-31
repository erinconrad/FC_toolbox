function [p z all_mu] = test_pt_circular_means(thing,polar_edges,hours_mins)

% Get circular mean for each patient
npts = size(thing,1);
nbins = size(thing,2);

% initialize rebin
mean_bins = zeros(size(thing));

all_mu = nan(npts,1);
for i = 1:npts
    curr = thing(i,:);
    all_mu(i) = wrapTo2Pi(circ_mean(polar_edges(1:nbins),curr'));
    
    for j = 1:nbins
        if all_mu(i) >= polar_edges(j) && all_mu(i) < polar_edges(j+1)
            mean_bins(i,j) = 1;
        end
    end
end

%assert(isequal(sum(mean_bins,2),ones(npts,1)));

sum_bins = nansum(mean_bins,1);
all_mu(isnan(all_mu)) = [];
[p z] = circ_rtest(all_mu);



if 0
figure
nexttile
circ_plot(all_mu,'hist',[],length(polar_edges),true,true,'linewidth',2,'color','r')

skip = 6;
nexttile
polarhistogram('BinEdges',polar_edges,'BinCounts',sum_bins,...
    'displayStyle','stairs','linewidth',2)
set(gca,'ThetaDir','clockwise');
set(gca,'ThetaZeroLocation','top');
set(gca,'rticklabels',[])
thetaticks(polar_edges(1:skip:nbins)*360/(2*pi))
thetaticklabels(hours_mins(1:skip:nbins+1))
set(gca,'fontsize',15)
title('Normalized spike rate by time of day')

end
%}